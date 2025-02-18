#%% 
import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
from pathlib import Path
from scipy.spatial import cKDTree
import argparse
import matplotlib.pyplot as plt
from dotenv import load_dotenv
from nibabel.freesurfer.io import read_geometry
import trimesh
from nilearn import plotting as niplot
from matplotlib.colors import LinearSegmentedColormap
import plotly.graph_objects as go

#%% 
class IEEGRecon:
    def __init__(self, pre_implant_mri, post_implant_ct, ct_electrodes, output_dir, env_path=None, freesurfer_dir=None):
        """
        Initialize IEEGRecon with required paths
        
        Args:
            pre_implant_mri (str): Path to pre-implant MRI
            post_implant_ct (str): Path to post-implant CT
            ct_electrodes (str): Path to electrode coordinates CSV
            output_dir (str): Output directory path (required)
            env_path (str, optional): Path to .env file
            freesurfer_dir (str, optional): Path to FreeSurfer subjects directory. If provided, overrides SUBJECTS_DIR from env
        """
        # Set main properties
        self.preImplantMRI = pre_implant_mri
        self.postImplantCT = post_implant_ct
        self.postImplantCT_electrodes = ct_electrodes
        self.output = Path(output_dir)  # Convert to Path object
        self.root_dir = Path(__file__).parent.parent
        
        # Setup environment variables and paths
        self._setup_environment(env_path, freesurfer_dir)

    def _setup_environment(self, env_path=None, freesurfer_dir=None):
        """
        Setup environment variables and paths from .env file
        
        Args:
            env_path (str, optional): Path to .env file
            freesurfer_dir (str, optional): Path to FreeSurfer subjects directory
        """
        if env_path is None:
            env_path = Path(__file__).parent / '.env'
        
        try:
            # Load environment variables from .env file
            load_dotenv(env_path)
            
            # Set paths from environment variables
            self.fslLoc = os.getenv('FSL_DIR')
            self.itksnap = os.getenv('ITKSNAP_DIR')
            self.freeSurfer = os.getenv('FREESURFER_HOME')
            # Allow freesurfer_dir parameter to override environment variable
            self.freeSurferDir = freesurfer_dir if freesurfer_dir is not None else os.getenv('SUBJECTS_DIR')
            
            # Ensure required environment variables are set
            if not all([self.fslLoc, self.itksnap, self.freeSurfer]):
                raise ValueError("Missing required environment variables")
            
            # Set FSL output type
            os.environ['FSLOUTPUTTYPE'] = os.getenv('FSL_OUTPUT_TYPE', 'NIFTI_GZ')
            
            # Run FreeSurfer setup if available
            freesurfer_setup = Path(self.freeSurfer) / 'SetUpFreeSurfer.sh'
            if freesurfer_setup.exists():
                subprocess.run(['sh', str(freesurfer_setup)], check=True)
                
        except Exception as e:
            print(f"Warning: Environment setup error - {str(e)}")
            print("Using default None values for paths")
            self.fslLoc = self.itksnap = self.freeSurfer = self.freeSurferDir = None

    def module1(self):
        """
        Module1: exports electrode coordinates of post implant CT in voxel and
        native space. Outputs of this module goes in output:ieeg_recon/module1 folder
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module1'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Export electrode coordinates in CT space in mm and vox
        elecCTvox = pd.read_csv(self.postImplantCT_electrodes,  sep=r'\s+', header=None)
            
        # Write electrode names
        with open(output_dir / 'electrode_names.txt', 'w') as f:
            f.write('\n'.join(elecCTvox.iloc[:, 0].tolist()))
        
        # Write electrode coordinates
        np.savetxt(
            output_dir / 'electrodes_inCTvox.txt',
            elecCTvox.iloc[:, 1:4].values.astype(int),  # Convert to integers
            delimiter=' ',
            fmt='%d'  # Use integer format
        )

        # Load CT data
        ct_img = nib.load(self.postImplantCT)
        ct_data = ct_img.get_fdata()
        ct_affine = ct_img.affine

        # Convert electrode coordinates from voxel to world space
        elec_vox = elecCTvox.iloc[:, 1:4].values
        elec_homog = np.hstack((elec_vox, np.ones((elec_vox.shape[0], 1))))  # Add homogeneous coordinate
        elecCTmm = np.dot(ct_affine, elec_homog.T).T[:, :3]  # Transform and remove homogeneous coordinate

        # Save world coordinates
        np.savetxt(
            output_dir / 'electrodes_inCTmm.txt',
            elecCTmm,
            delimiter=' ',
            fmt='%.2f'  # Use float format
        )

    def module2(self, reg_type, skip_existing=False):
        """
        Module2: Outputs go in output:ieeg_recon/module2 folder
        
        Args:
            reg_type (str): Registration type - 'gc', 'g', or 'gc_noCTthereshold'
            skip_existing (bool): If True, skip processing if output files exist
        
        Returns:
            dict: Paths to output files
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module2'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Define output file locations
        file_locations = {
            'ct_to_mri': output_dir / 'ct_to_mri.nii.gz',
            'electrodes_inMRI': output_dir / 'electrodes_inMRI.nii.gz',
            'electrodes_inMRI_freesurferLUT': output_dir / 'electrodes_inMRI_freesurferLUT.txt',
            'electrodes_inMRImm': output_dir / 'electrodes_inMRImm.txt',
            'electrodes_inMRIvox': output_dir / 'electrodes_inMRIvox.txt'
        }

        # Check if files exist and skip if requested
        if skip_existing and all(path.exists() for path in file_locations.values()):
            return file_locations

        # Remove negative values from CT image
        cmd = ['fslmaths', 
               self.postImplantCT, 
               '-thr', '0', 
               str(output_dir / 'ct_thresholded.nii.gz')]
        subprocess.run(cmd, check=True)

        # Handle different registration types
        if reg_type == 'gc':
            self._run_greedy_centered_registration(output_dir)
        elif reg_type == 'g':
            self._run_flirt_greedy_registration(output_dir)
        elif reg_type == 'gc_noCTthereshold':
            self._run_greedy_centered_no_threshold(output_dir)
        else:
            raise ValueError("Registration type must be 'gc', 'g', or 'gc_noCTthereshold'")

        # Apply registration to electrode coordinates
        self._transform_electrode_coordinates(output_dir)
        
        # Create electrode spheres
        self._create_electrode_spheres(output_dir)

        return file_locations

    def _run_greedy_centered_registration(self, output_dir):
        """Run greedy registration with image centering"""
        print('Running greedy registration with image centering')
        
        # Run greedy registration
        subprocess.run([
            f"{self.itksnap}/greedy",
            "-d", "3",
            "-i", self.preImplantMRI,
            str(output_dir / 'ct_thresholded.nii.gz'),
            "-o", str(output_dir / 'ct_to_mri.mat'),
            "-a", "-dof", "6",
            "-m", "NMI",
            "-ia-image-centers",
            "-n", "100x100x0x0",
            "-jitter", "0.9",
            "-search", "1000", "10", "20"
        ], stdout=open(output_dir / 'greedy.log', 'w'), check=True)

        # Convert transform to FSL format
        subprocess.run([
            f"{self.itksnap}/c3d_affine_tool",
            "-ref", self.preImplantMRI,
            "-src", str(output_dir / 'ct_thresholded.nii.gz'),
            str(output_dir / 'ct_to_mri.mat'),
            "-ras2fsl",
            "-o", str(output_dir / 'ct_to_mri_xform.txt')
        ], check=True)

        # Remove temporary mat file
        (output_dir / 'ct_to_mri.mat').unlink()

        # Apply transform
        cmd = ['flirt',
               '-in', str(output_dir / 'ct_thresholded.nii.gz'),
               '-ref', self.preImplantMRI,
               '-init', str(output_dir / 'ct_to_mri_xform.txt'),
               '-out', str(output_dir / 'ct_to_mri.nii.gz'),
               '-applyxfm']
        subprocess.run(cmd, check=True)

    def _run_greedy_centered_no_threshold(self, output_dir):
        """Run greedy registration with image centering without CT thresholding"""
        print('Running greedy registration with image centering')
        
        # Run greedy registration
        subprocess.run([
            f"{self.itksnap}/greedy",
            "-d", "3",
            "-i", self.preImplantMRI,
            self.postImplantCT,  # Use original CT instead of thresholded
            "-o", str(output_dir / 'ct_to_mri.mat'),
            "-a", "-dof", "6",
            "-m", "NMI",
            "-ia-image-centers",
            "-n", "100x50x0x0"
        ], stdout=open(output_dir / 'greedy.log', 'w'), check=True)

        # Convert transform to FSL format
        subprocess.run([
            f"{self.itksnap}/c3d_affine_tool",
            "-ref", self.preImplantMRI,
            "-src", self.postImplantCT,
            str(output_dir / 'ct_to_mri.mat'),
            "-ras2fsl",
            "-o", str(output_dir / 'ct_to_mri_xform.txt')
        ], check=True)

        # Remove temporary mat file
        (output_dir / 'ct_to_mri.mat').unlink()

        # Apply transform
        cmd = ['flirt',
               '-in', self.postImplantCT,
               '-ref', self.preImplantMRI,
               '-init', str(output_dir / 'ct_to_mri_xform.txt'),
               '-out', str(output_dir / 'ct_to_mri.nii.gz'),
               '-applyxfm']
        subprocess.run(cmd, check=True)

        # Threshold the registered CT image
        cmd = ['fslmaths',
               str(output_dir / 'ct_to_mri.nii.gz'),
               '-thr', '0',
               str(output_dir / 'ct_to_mri.nii.gz')]
        subprocess.run(cmd, check=True)

    def _transform_electrode_coordinates(self, output_dir):
        """Apply registration transform to electrode coordinates"""
        # Transform mm coordinates
        subprocess.run(["img2imgcoord",
            "-src", str(output_dir / 'ct_thresholded.nii.gz'),
            "-dest", str(output_dir / 'ct_to_mri.nii.gz'),
            "-xfm", str(output_dir / 'ct_to_mri_xform.txt'),
            "-mm", str(Path(self.output) / 'ieeg_recon/module1/electrodes_inCTmm.txt')
        ], stdout=open(output_dir / 'electrodes_inMRImm.txt', 'w'), check=True)

        # Transform voxel coordinates
        subprocess.run(["img2imgcoord",
            "-src", str(output_dir / 'ct_thresholded.nii.gz'),
            "-dest", str(output_dir / 'ct_to_mri.nii.gz'),
            "-xfm", str(output_dir / 'ct_to_mri_xform.txt'),
            "-vox", str(Path(self.output) / 'ieeg_recon/module1/electrodes_inCTvox.txt')
        ], stdout=open(output_dir / 'electrodes_inMRIvox.txt', 'w'), check=True)

    def _create_electrode_spheres(self, output_dir):
        """Create spheres for electrodes in registered space"""
        # Load registered CT data
        ct_img = nib.load(output_dir / 'ct_to_mri.nii.gz')
        ct_data = ct_img.get_fdata()
        ct_affine = ct_img.affine

        # Create blank image
        blank_data = np.zeros_like(ct_data)
        vox_coords = np.array(np.where(blank_data == 0)).T
        
        # Convert to world coordinates
        vox_homog = np.hstack((vox_coords, np.ones((vox_coords.shape[0], 1))))
        world_coords = np.dot(ct_affine, vox_homog.T).T[:, :3]

        # Load electrode coordinates
        electrodes_mm = np.loadtxt(output_dir / 'electrodes_inMRImm.txt', skiprows=1)
        electrode_names = np.loadtxt(
            Path(self.output) / 'ieeg_recon/module1/electrode_names.txt',
            dtype=str
        )

        # Create FreeSurfer LUT
        n_electrodes = len(electrode_names)
        lut_data = {
            'index': np.arange(1, n_electrodes + 1),
            'names': electrode_names,
            'R': np.full(n_electrodes, 90),
            'G': np.full(n_electrodes, 150),
            'B': np.full(n_electrodes, 60),
            'alpha': np.zeros(n_electrodes)
        }
        pd.DataFrame(lut_data).to_csv(
            output_dir / 'electrodes_inMRI_freesurferLUT.txt',
            sep=' ',
            header=False,
            index=False
        )

        # Find points within 2mm of each electrode
        tree = cKDTree(world_coords)
        dist, idx = tree.query(electrodes_mm, k=1)
        
        # Create electrode map
        electrode_data = blank_data.copy()
        mask = dist <= 2
        for i, (valid, coord) in enumerate(zip(mask, vox_coords[idx]), 1):
            if valid:
                electrode_data[tuple(coord)] = i

        # Save electrode map
        nib.save(
            nib.Nifti1Image(electrode_data, ct_affine),
            output_dir / 'electrodes_inMRI.nii.gz'
        )

    def module2_QualityAssurance(self, file_locations, imageviewer):
        """
        Generate quality assurance visualizations for module2 results
        
        Args:
            file_locations (dict): Dictionary containing paths to module2 output files
            imageviewer (str): Type of viewer to use - 'freeview_snapshot', 'freeview', 'niplot', or 'itksnap'
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module2'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Only launch viewer if not set to 'none'
        if imageviewer.lower() != 'none':
            if imageviewer == 'freeview':
                # Open interactive freeview session
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "sagittal"
                ], check=True)
            elif imageviewer == 'freeview_snapshot':
                # Generate sagittal view
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "sagittal",
                    "-ss", str(output_dir / "QA_registation_sagittal.png")
                ], check=True)

                # Generate coronal view
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "coronal",
                    "-ss", str(output_dir / "QA_registation_coronal.png")
                ], check=True)

                # Generate axial view
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "axial",
                    "-ss", str(output_dir / "QA_registation_axial.png")
                ], check=True)

                # Generate 3D view
                subprocess.run([
                    "freeview",
                    "-v", f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}:isosurface=on",
                    "-viewport", "3d", "-view", "anterior",
                    "-ss", str(output_dir / "QA_registation_3D.png")
                ], check=True)
            elif imageviewer == 'itksnap':
                # Open interactive ITK-SNAP session
                subprocess.run([
                    f"{self.itksnap}/itksnap",
                    "-g", self.preImplantMRI,
                    "-o", file_locations['ct_to_mri']
                ], check=True)
            elif imageviewer == 'niplot':
                # Create custom colormap that is transparent for zeros and scales from yellow to red
                colors = [(0, 0, 0, 0),          # transparent
                         (1, 1, 0, 1),           # yellow
                         (1, 0.5, 0, 1),         # orange
                         (1, 0, 0, 1)]           # red
                n_bins = 256  # Number of gradients
                custom_cmap = LinearSegmentedColormap.from_list('custom_hot', colors, N=n_bins)
                
                # Create the plot
                display = niplot.plot_roi(
                    file_locations['ct_to_mri'], 
                    bg_img=self.preImplantMRI, 
                    cmap=custom_cmap,
                    title='Registration', 
                    display_mode="mosaic"
                )
                
                # Save the plot
                plt.savefig(str(output_dir / "QA_registration_niplot.png"), 
                           dpi=300, 
                           bbox_inches='tight')
                
                # Show the plot and wait for it to be closed
                plt.show(block=False)
            else:
                raise ValueError(f"Unknown imageviewer option: {imageviewer}")

    def module3(self, atlas, atlas_lut, diameter=2.5, skip_existing=False):
        """
        Module3: Map electrodes to brain regions using provided atlas
        
        Args:
            atlas (str/Path): Path to atlas NIFTI file
            atlas_lut (str/Path): Path to lookup table CSV/txt file
            diameter (float): Maximum distance in mm for electrode-to-ROI mapping (default: 2.5)
            skip_existing (bool): If True, skip processing if output files exist
        
        Returns:
            str: Path to output electrodes2ROI CSV file
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon/module3'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define output file locations
        file_locations = {
            'electrodes2ROI': output_dir / 'electrodes2ROI.csv'
        }

        # Check if files exist and skip if requested
        if skip_existing and all(path.exists() for path in file_locations.values()):
            return str(file_locations['electrodes2ROI'])

        # Load electrode coordinates and names
        electrodes_mm = np.loadtxt(
            Path(self.output) / 'ieeg_recon/module2/electrodes_inMRImm.txt',
            skiprows=1
        )
        
        electrodes_vox = np.loadtxt(
            Path(self.output) / 'ieeg_recon/module2/electrodes_inMRIvox.txt',
            skiprows=1
        )
        
        labels = np.loadtxt(
            Path(self.output) / 'ieeg_recon/module1/electrode_names.txt',
            dtype=str
        )

        # Load atlas and lookup table
        atlas_img = nib.load(atlas)
        atlas_data = atlas_img.get_fdata()
        lut = pd.read_csv(atlas_lut, sep=None, engine='python')

        # Get vox2ras-tkr transform from atlas header
        # vox2ras = atlas_img.header.get_vox2ras()
        vox2ras_tkr = atlas_img.header.get_vox2ras_tkr()
        
        # Transform electrode coordinates to surface space
        electrodes_homog = np.hstack((electrodes_vox, np.ones((electrodes_vox.shape[0], 1))))
        electrodes_surfmm = np.round(np.dot(vox2ras_tkr, electrodes_homog.T).T[:, :3], decimals=4)

        # Get atlas ROI coordinates
        atlas_voxels = []
        for _, row in lut.iterrows():
            vox = np.array(np.where(atlas_data == row['roiNum'])).T
            if len(vox) > 0:
                atlas_voxels.append(np.column_stack([vox, np.full(len(vox), row['roiNum'])]))
        
        atlas_voxels = np.vstack(atlas_voxels)

        # Convert atlas voxels to mm space
        vox_homog = np.hstack((atlas_voxels[:, :3], np.ones((len(atlas_voxels), 1))))
        cord_mm = np.dot(atlas_img.affine, vox_homog.T).T[:, :3]
        
        # Find nearest ROI for each electrode
        tree = cKDTree(cord_mm)
        dist_mm, idx = tree.query(electrodes_mm, k=1)
        
        # Get ROI numbers for each electrode
        implant2roiNum = atlas_voxels[idx, 3]
        
        # Map ROI numbers to names
        implant2roi = pd.Series(implant2roiNum).map(
            lut.set_index('roiNum')['roi'].to_dict()
        ).fillna('')
        
        # Mark contacts beyond diameter as white matter or outside brain
        implant2roi[dist_mm > diameter] = ''
        implant2roiNum = pd.Series(implant2roiNum)
        implant2roiNum[dist_mm > diameter] = pd.NA

        # Load FreeSurfer surfaces
        lh_pial_verts, lh_pial_faces = read_geometry(Path(self.freeSurferDir) / 'surf/lh.pial')
        lh_white_verts, lh_white_faces = read_geometry(Path(self.freeSurferDir) / 'surf/lh.white')
        rh_pial_verts, rh_pial_faces = read_geometry(Path(self.freeSurferDir) / 'surf/rh.pial')
        rh_white_verts, rh_white_faces = read_geometry(Path(self.freeSurferDir) / 'surf/rh.white')

        # Create surface meshes
        lh_pial_mesh = trimesh.Trimesh(vertices=lh_pial_verts, faces=lh_pial_faces)
        rh_pial_mesh = trimesh.Trimesh(vertices=rh_pial_verts, faces=rh_pial_faces)
        lh_white_mesh = trimesh.Trimesh(vertices=lh_white_verts, faces=lh_white_faces)
        rh_white_mesh = trimesh.Trimesh(vertices=rh_white_verts, faces=rh_white_faces)

        # Find contacts outside brain
        outside_mask = ~(
            lh_pial_mesh.contains(electrodes_surfmm[dist_mm > diameter]) | 
            rh_pial_mesh.contains(electrodes_surfmm[dist_mm > diameter])
        )
        outside_indices = np.where(dist_mm > diameter)[0][outside_mask]
        implant2roi.iloc[outside_indices] = 'outside-brain'

        # Find white matter contacts
        wm_indices = np.where(dist_mm > diameter)[0][~outside_mask]
        implant2roi.iloc[wm_indices] = 'white-matter'

        # Verify white matter contacts - check if contacts labeled as white matter
        # are actually inside the white matter surface
        wm_contacts = electrodes_surfmm[wm_indices]
        in_left = lh_white_mesh.contains(wm_contacts)
        in_right = rh_white_mesh.contains(wm_contacts)
        in_white = in_left | in_right

        if np.sum(in_white) == len(wm_indices):
            print('white-matter contacts correctly assigned')
        else:
            print(f"check white-matter contacts in: {self.output}")

        # Create output dataframe
        electrodes2ROI = pd.DataFrame({
            'labels': labels,
            'mm_x': electrodes_mm[:, 0],
            'mm_y': electrodes_mm[:, 1],
            'mm_z': electrodes_mm[:, 2],
            'surfmm_x': electrodes_surfmm[:, 0],
            'surfmm_y': electrodes_surfmm[:, 1],
            'surfmm_z': electrodes_surfmm[:, 2],
            'vox_x': electrodes_vox[:, 0],
            'vox_y': electrodes_vox[:, 1],
            'vox_z': electrodes_vox[:, 2],
            'roi': implant2roi,
            'roiNum': implant2roiNum
        })

        # Save output
        electrodes2ROI.to_csv(file_locations['electrodes2ROI'], index=False)
        return str(file_locations['electrodes2ROI'])
    
    def module3_QualityAssurance(self, recon_file):
        """
        Generate quality assurance visualizations for module3 results
        
        Args:
            file_locations (str): Path to module3 output file
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module3'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load electrode data
        ieeg = pd.read_csv(recon_file)

        # Load FreeSurfer surfaces
        lh_pial_verts, lh_pial_faces = read_geometry(Path(self.freeSurferDir) / 'surf/lh.pial')
        rh_pial_verts, rh_pial_faces = read_geometry(Path(self.freeSurferDir) / 'surf/rh.pial')

        # Create interactive 3D visualization using plotly
        fig = go.Figure()

        # Update mesh settings for left hemisphere
        fig.add_trace(go.Mesh3d(
            x=lh_pial_verts[:, 0], y=lh_pial_verts[:, 1], z=lh_pial_verts[:, 2],
            i=lh_pial_faces[:, 0], j=lh_pial_faces[:, 1], k=lh_pial_faces[:, 2],
            color='#808080',  # Medium grey color
            opacity=1.0,
            flatshading=False,
            lighting=dict(
                ambient=0.3,     # Reduced ambient light for more contrast
                diffuse=0.8,     # Strong diffuse for material look
                specular=1.0,    # Maximum specular for shininess
                roughness=0.1,   # Low roughness for smoother, shinier surface
                fresnel=0.9      # High fresnel for metallic appearance
            ),
            lightposition=dict(x=100, y=200, z=150),  # Angled light for better highlights
            name='Left Hemisphere',
            visible=True,
            showscale=False
        ))

        # Update mesh settings for right hemisphere (same settings as left)
        fig.add_trace(go.Mesh3d(
            x=rh_pial_verts[:, 0], y=rh_pial_verts[:, 1], z=rh_pial_verts[:, 2],
            i=rh_pial_faces[:, 0], j=rh_pial_faces[:, 1], k=rh_pial_faces[:, 2],
            color='#808080',  # Medium grey color
            opacity=1.0,
            flatshading=False,
            lighting=dict(
                ambient=0.3,     # Reduced ambient light for more contrast
                diffuse=0.8,     # Strong diffuse for material look
                specular=1.0,    # Maximum specular for shininess
                roughness=0.1,   # Low roughness for smoother, shinier surface
                fresnel=0.9      # High fresnel for metallic appearance
            ),
            lightposition=dict(x=100, y=200, z=150),  # Angled light for better highlights
            name='Right Hemisphere',
            visible=True,
            showscale=False
        ))

        # Update electrode points
        fig.add_trace(go.Scatter3d(
            x=ieeg['surfmm_x'],
            y=ieeg['surfmm_y'],
            z=ieeg['surfmm_z'],
            mode='markers',
            marker=dict(
                size=4,
                color='red',
                opacity=1.0,
                symbol='circle'
            ),
            text=ieeg['labels'],
            hovertemplate=(
                "<b>%{text}</b><br>" +
                "ROI: %{customdata[0]}<br>" +
                "X: %{customdata[1]:.2f}<br>" +
                "Y: %{customdata[2]:.2f}<br>" +
                "Z: %{customdata[3]:.2f}<br>" +
                "<extra></extra>"
            ),
            customdata=np.column_stack((
                ieeg['roi'],
                ieeg['surfmm_x'],
                ieeg['surfmm_y'],
                ieeg['surfmm_z']
            )),
            name='Electrodes',
            showlegend=True
        ))

        # Create buttons for toggling hemisphere visibility and electrode labels
        updatemenus = [
            # Hemisphere visibility buttons - Top
            dict(
                type="buttons",
                showactive=True,
                buttons=[
                    dict(label="Show All",
                         method="update",
                         args=[{"visible": [True, True, True]}]),
                    dict(label="Left Only",
                         method="update",
                         args=[{"visible": [True, False, True]}]),
                    dict(label="Right Only",
                         method="update",
                         args=[{"visible": [False, True, True]}]),
                    dict(label="Electrodes Only",
                         method="update",
                         args=[{"visible": [False, False, True]}]),
                ],
                direction="down",
                pad={"r": 10, "t": 10},
                x=0.02,  # Move closer to left edge
                xanchor="left",
                y=0.9,   # Position at top
                yanchor="top"
            ),
            # Electrode size buttons - Middle
            dict(
                type="buttons",
                showactive=True,
                buttons=[
                    dict(label="Small Electrodes",
                         method="restyle",
                         args=[{"marker.size": 5}]),
                    dict(label="Medium Electrodes",
                         method="restyle",
                         args=[{"marker.size": 8}]),
                    dict(label="Large Electrodes",
                         method="restyle",
                         args=[{"marker.size": 10}]),
                ],
                direction="down",
                pad={"r": 10, "t": 10},
                x=0.02,  # Move closer to left edge
                xanchor="left",
                y=0.5,   # Position in middle
                yanchor="middle"
            ),
            # Label visibility toggle - Bottom
            dict(
                type="buttons",
                showactive=True,
                buttons=[
                    dict(label="Show Labels",
                         method="restyle",
                         args=[{"mode": "markers+text"}]),
                    dict(label="Hide Labels",
                         method="restyle",
                         args=[{"mode": "markers"}]),
                ],
                direction="down",
                pad={"r": 10, "t": 10},
                x=0.02,  # Move closer to left edge
                xanchor="left",
                y=0.1,   # Position at bottom
                yanchor="bottom"
            )
        ]

        # Update layout with brighter settings
        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False, showgrid=False, showbackground=False),
                yaxis=dict(visible=False, showgrid=False, showbackground=False),
                zaxis=dict(visible=False, showgrid=False, showbackground=False),
                aspectmode="data",
                camera=dict(
                    up=dict(x=0, y=0, z=1),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=1.5, y=1.5, z=1.5)
                ),
                bgcolor='white'
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(r=0, l=0, b=0, t=0),
            showlegend=False,
            updatemenus=updatemenus
        )

        # Save the plot as HTML for interactive viewing
        fig.write_html(str(output_dir / 'electrode_visualization.html'))

        # Save as a static image at 300 DPI
        fig.write_image(str(output_dir / 'electrode_visualization.png'), scale=3)

#%%
def run_pipeline(pre_implant_mri, 
                 post_implant_ct, 
                 ct_electrodes, 
                 output_dir, 
                 env_path=None, 
                 freesurfer_dir=None,
                 modules=['1', '2', '3'], 
                 skip_existing=False, 
                 reg_type='gc_noCTthereshold', 
                 qa_viewer='niplot'):
    """
    Run the iEEG reconstruction pipeline
    
    Args:
        pre_implant_mri (str): Path to pre-implant MRI
        post_implant_ct (str): Path to post-implant CT
        ct_electrodes (str): Path to electrode coordinates CSV
        output_dir (str/Path): Output directory path
        env_path (str/Path, optional): Path to .env file
        modules (list): List of modules to run ['1', '2']
        skip_existing (bool): Skip processing if output files exist
        reg_type (str): Registration type ('gc', 'g', 'gc_noCTthereshold')
        qa_viewer (str): Quality assurance viewer type
    
    Returns:
        dict: Paths to output files (if module 2 was run)
    """
    # Initialize reconstruction object
    recon = IEEGRecon(
        pre_implant_mri=pre_implant_mri,
        post_implant_ct=post_implant_ct,
        ct_electrodes=ct_electrodes,
        output_dir=output_dir,
        env_path=env_path,
        freesurfer_dir=freesurfer_dir
    )
    
    # Run selected modules
    file_locations = None
    
    if '1' in modules:
        print("Running Module 1...")
        recon.module1()
    
    if '2' in modules:
        print("Running Module 2...")
        file_locations = recon.module2(reg_type, skip_existing=skip_existing)
        
        print("Output files:")
        for name, path in file_locations.items():
            print(f"{name}: {path}")
        
        recon.module2_QualityAssurance(file_locations, qa_viewer)

    if '3' in modules:
        print("Running Module 3...")
        atlas = freesurfer_dir / 'mri' / 'aparc+aseg.mgz'
        atlas_lut = project_path / 'doc' / 'atlasLUT' / 'desikanKilliany.csv'
        recon_file = recon.module3(atlas, atlas_lut, diameter=2.5, skip_existing=skip_existing)
        recon.module3_QualityAssurance(recon_file)
    
    return file_locations

#%%
if __name__ == "__main__":
    # Example usage - replace these values with your actual file paths
    project_path = Path(__file__).parent.parent
    filepath_csv = project_path / 'test' / 'test.csv'
    
    # Read filepath CSV
    filepath = pd.read_csv(filepath_csv, index_col='record_id')
       
    # Set paths for the selected subject
    pre_implant_mri = filepath.loc['sub-RID0031']['t1w']
    post_implant_ct = filepath.loc['sub-RID0031']['ct']
    ct_electrodes = filepath.loc['sub-RID0031']['electrodes']
    output_dir = project_path / 'test' / 'output'
    freesurfer_dir = project_path / 'data' / 'sub-RID0031' / 'derivatives' / 'freesurfer'
   
    # Set config path (defaults to .env in same directory as script)
    env_path = project_path / '.env'
    
    print(f"Processing test subject")
    
    # Run pipeline with default settings
    run_pipeline(
        pre_implant_mri=pre_implant_mri,
        post_implant_ct=post_implant_ct,
        ct_electrodes=ct_electrodes,
        output_dir=output_dir,
        env_path=env_path,
        freesurfer_dir=freesurfer_dir,
        modules=['1', '2', '3'],
        skip_existing=True,
        reg_type='gc_noCTthereshold',  # Default registration type
        qa_viewer='niplot'  # Default viewer
    )
    
    print("Processing complete!")

# %%
