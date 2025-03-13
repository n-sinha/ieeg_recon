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
import plotly.express as px
from IPython import embed

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

        embed()

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
        ).astype(int)
        
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

    def module4(self, skip_existing=False):
        """
        Module4: Transform electrode coordinates to MNI305 and MNI152 spaces
        
        Args:
            skip_existing (bool): If True, skip processing if output files exist
        
        Returns:
            dict: Paths to output files
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module4'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Define output file locations
        file_locations = {
            'electrodes2ROI_mni305': output_dir / 'electrodes2ROI_mni305.csv',
            'electrodes2ROI_mni152': output_dir / 'electrodes2ROI_mni152.csv'
        }

        # Check if files exist and skip if requested
        if skip_existing and all(path.exists() for path in file_locations.values()):
            return file_locations

        # Define input paths
        recon_native = self.output / 'ieeg_recon' / 'module3' / 'electrodes2ROI.csv'
        mni305 = Path(self.freeSurfer) / 'subjects' / 'fsaverage'
        mni152 = Path(self.freeSurfer) / 'subjects' / 'cvs_avg35_inMNI152'
        talXFM_path = Path(self.freeSurferDir) / 'mri' / 'transforms' / 'talairach.xfm'
        t1mgz_path = Path(self.freeSurferDir) / 'mri' / 'T1.mgz'

        # Verify input files exist
        required_paths = {
            'Module 3 output': recon_native,
            'MNI 305 template': mni305,
            'MNI 152 template': mni152,
            'Talairach transform': talXFM_path,
            'T1 MGZ file': t1mgz_path
        }
        for name, path in required_paths.items():
            if not path.exists():
                raise FileNotFoundError(f"{name} not found: {path}")

        # Load electrode data
        electrodes2ROI = pd.read_csv(recon_native)

        # Get all transformations
        # ref: https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
        xform_tal = self._read_talairach_xfm(talXFM_path)
        t1mgz = nib.load(t1mgz_path)
        Norig = t1mgz.header.get_vox2ras()
        Torig = t1mgz.header.get_vox2ras_tkr()
        xform_mni305 = np.dot(xform_tal, np.dot(Norig, np.linalg.inv(Torig)))

        # Load template spaces
        mni305_t1mgz = nib.load(mni305 / 'mri' / 'T1.mgz')
        mni152_t1mgz = nib.load(mni152 / 'mri' / 'T1.mgz')
        xform_mni305_tk_ras = mni305_t1mgz.header.get_vox2ras_tkr()
        xform_mni152_tk_ras = mni152_t1mgz.header.get_vox2ras_tkr()

        # Define MNI152 transform matrix
        xform_mni152 = np.array([
            [0.9975, -0.0073, 0.0176, -0.0429],
            [0.0146, 1.0009, -0.0024, 1.5496],
            [-0.0130, -0.0093, 0.9971, 1.1840],
            [0.0, 0.0, 0.0, 1.0]
        ])

        # Transform coordinates to MNI305 space
        surfmm = electrodes2ROI.filter(['surfmm_x', 'surfmm_y', 'surfmm_z']).to_numpy()
        surfmm_homog = np.hstack((surfmm, np.ones((surfmm.shape[0], 1))))
        mni305_surfmm = np.round(np.dot(xform_mni305, surfmm_homog.T).T[:, :3], decimals=4)
        mni305_vox = np.dot(np.linalg.inv(xform_mni305_tk_ras),
                           np.hstack((mni305_surfmm, np.ones((mni305_surfmm.shape[0], 1)))).T).T[:, :3].astype(int)
        mni305_mm = nib.affines.apply_affine(mni305_t1mgz.affine, mni305_vox)

        # Transform coordinates to MNI152 space
        mni305_surfmm_homog = np.hstack((mni305_surfmm, np.ones((mni305_surfmm.shape[0], 1))))
        mni152_surfmm = np.round(np.dot(xform_mni152, mni305_surfmm_homog.T).T[:, :3], decimals=4)
        mni152_vox = np.dot(np.linalg.inv(xform_mni152_tk_ras), 
                           np.hstack((mni152_surfmm, np.ones((mni152_surfmm.shape[0], 1)))).T).T[:, :3].astype(int)
        mni152_mm = nib.affines.apply_affine(mni152_t1mgz.affine, mni152_vox)

        # Create and save MNI305 coordinates DataFrame
        electrodes2ROI_mni305 = pd.DataFrame({
            'labels': electrodes2ROI['labels'],
            'mni305_mm_x': mni305_mm[:, 0],
            'mni305_mm_y': mni305_mm[:, 1],
            'mni305_mm_z': mni305_mm[:, 2],
            'mni305_surfmm_x': mni305_surfmm[:, 0],
            'mni305_surfmm_y': mni305_surfmm[:, 1],
            'mni305_surfmm_z': mni305_surfmm[:, 2],
            'mni305_vox_x': mni305_vox[:, 0],
            'mni305_vox_y': mni305_vox[:, 1],
            'mni305_vox_z': mni305_vox[:, 2],
            'roi': electrodes2ROI['roi'],
            'roiNum': electrodes2ROI['roiNum']
        })
        electrodes2ROI_mni305.to_csv(file_locations['electrodes2ROI_mni305'], index=False)

        # Create and save MNI152 coordinates DataFrame
        electrodes2ROI_mni152 = pd.DataFrame({
            'labels': electrodes2ROI['labels'],
            'mni152_mm_x': mni152_mm[:, 0],
            'mni152_mm_y': mni152_mm[:, 1],
            'mni152_mm_z': mni152_mm[:, 2],
            'mni152_surfmm_x': mni152_surfmm[:, 0],
            'mni152_surfmm_y': mni152_surfmm[:, 1],
            'mni152_surfmm_z': mni152_surfmm[:, 2],
            'mni152_vox_x': mni152_vox[:, 0],
            'mni152_vox_y': mni152_vox[:, 1],
            'mni152_vox_z': mni152_vox[:, 2],
            'roi': electrodes2ROI['roi'],
            'roiNum': electrodes2ROI['roiNum']
        })
        electrodes2ROI_mni152.to_csv(file_locations['electrodes2ROI_mni152'], index=False)

        return file_locations

    def _read_talairach_xfm(self, fname):
        """Read the transformation matrix from a FreeSurfer .xfm file.
        
        Args:
            fname (str/Path): Path to FreeSurfer .xfm file
            
        Returns:
            numpy.ndarray: 4x4 transformation matrix
        """
        # Skip header lines until we find 'Linear_Transform'
        xfm = []
        with open(fname) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'Linear_Transform' in line:
                    # Read the next 3 lines as the transformation matrix
                    for j in range(3):
                        numbers = [float(x) for x in lines[i + 1 + j].strip('\n;').split()]
                        xfm.append(numbers)
        
        # Add the last row [0, 0, 0, 1] to make it a 4x4 matrix
        xfm.append([0.0, 0.0, 0.0, 1.0])
        return np.array(xfm)
    
    def module4_snap_to_atlas(self, 
                              standard_space='mni152', 
                              atlas='aparc+aseg.mgz', 
                              atlas_lut='desikanKilliany.csv', 
                              diameter=2.5) -> Path:
        """
        Module4: Snap electrodes to atlas ROIs in standard space
        
        Args:
            standard_space (str): Standard space to use - must be either 'mni305' or 'mni152'
            atlas (str/Path): Path to atlas NIFTI file
            atlas_lut (str/Path): Path to lookup table CSV/txt file
            diameter (float): Maximum distance in mm for electrode-to-ROI mapping (default: 2.5)
            
        Raises:
            ValueError: If standard_space is not 'mni305' or 'mni152'
        """
        # Validate standard space input
        if standard_space.lower() not in ['mni305', 'mni152']:
            raise ValueError("standard_space must be either 'mni305' or 'mni152'")
        
        if standard_space == 'mni305':
            standard_file = self.output / 'ieeg_recon' / 'module4' / 'electrodes2ROI_mni305.csv'
            atlas_path = Path(self.freeSurfer) / 'subjects' / 'fsaverage'
        elif standard_space == 'mni152':
            standard_file = self.output / 'ieeg_recon' / 'module4' / 'electrodes2ROI_mni152.csv'
            atlas_path = Path(self.freeSurfer) / 'subjects' / 'cvs_avg35_inMNI152'

        # Load template spaces
        atlas_img = nib.load(atlas_path / 'mri' / 'aparc+aseg.mgz')
        atlas_data = atlas_img.get_fdata()
        lut = pd.read_csv(atlas_lut, sep=None, engine='python')
        recon_standard = pd.read_csv(standard_file)
        channels_to_snap = recon_standard[~recon_standard['roi'].isin(['white-matter', 'outside-brain'])]
        channels_to_skip = recon_standard[recon_standard['roi'].isin(['white-matter', 'outside-brain'])]

        # Get atlas ROI coordinates
        atlas_voxels = []
        for _, row in lut.iterrows():
            vox = np.array(np.where(atlas_data == row['roiNum'])).T
            if len(vox) > 0:
                atlas_voxels.append(np.column_stack([vox, np.full(len(vox), row['roiNum'])]))
        
        atlas_voxels = np.vstack(atlas_voxels)

        # Convert atlas voxels to mm space
        vox_homog = np.hstack((atlas_voxels[:, :3], np.ones((len(atlas_voxels), 1))))
        atlas_mm = np.dot(atlas_img.affine, vox_homog.T).T[:, :3]

        # Find nearest ROI for each electrode
        tree = cKDTree(atlas_mm)
        dist_mm, idx = tree.query(channels_to_snap.filter([f'{standard_space}_mm_x', 
                                                           f'{standard_space}_mm_y', 
                                                           f'{standard_space}_mm_z']), k=1)
        
         # Get ROI numbers for each electrode
        implant2roiNum = atlas_voxels[idx, 3].astype(int)
        # Create mask for channels to snap
        mask = channels_to_snap['roiNum'].astype(int) != implant2roiNum
        channels_to_snap = channels_to_snap[mask]

        for rowidx, row in channels_to_snap.iterrows():
            target_roi = atlas_mm[np.where(atlas_voxels[:,3] == row['roiNum'])]
            target_roi_tree = cKDTree(target_roi)
            dist_mm, idx = target_roi_tree.query(np.array([row[f'{standard_space}_mm_x'], 
                                               row[f'{standard_space}_mm_y'], 
                                               row[f'{standard_space}_mm_z']]), k=1)
            
            channels_to_snap.loc[rowidx, f'{standard_space}_mm_x'] = target_roi[idx, 0]
            channels_to_snap.loc[rowidx, f'{standard_space}_mm_y'] = target_roi[idx, 1]
            channels_to_snap.loc[rowidx, f'{standard_space}_mm_z'] = target_roi[idx, 2]

        # convert mm to voxels
        mm_homog = np.hstack((channels_to_snap.filter([f'{standard_space}_mm_x', 
                                                      f'{standard_space}_mm_y', 
                                                      f'{standard_space}_mm_z']), 
                              np.ones((len(channels_to_snap), 1))))
        
        voxels = np.dot(np.linalg.inv(atlas_img.affine), mm_homog.T).T[:, :3].astype(int)
        channels_to_snap[f'{standard_space}_vox_x'] = voxels[:, 0]
        channels_to_snap[f'{standard_space}_vox_y'] = voxels[:, 1]
        channels_to_snap[f'{standard_space}_vox_z'] = voxels[:, 2]

        xform_vox2ras = atlas_img.header.get_vox2ras()
        voxels_homog = np.hstack((channels_to_snap.filter([f'{standard_space}_vox_x', 
                                                          f'{standard_space}_vox_y', 
                                                          f'{standard_space}_vox_z']).to_numpy(), 
                                  np.ones((len(channels_to_snap), 1))))
        surf_mm = np.dot(xform_vox2ras, voxels_homog.T).T[:, :3]
        channels_to_snap[f'{standard_space}_mm_x'] = surf_mm[:, 0]
        channels_to_snap[f'{standard_space}_mm_y'] = surf_mm[:, 1]
        channels_to_snap[f'{standard_space}_mm_z'] = surf_mm[:, 2]

        # replace recon_standard with channels_to_snap where index matches
        recon_standard.loc[channels_to_snap.index, :] = channels_to_snap

        # Create new filename with 'corrected' in it
        standard_file_path = Path(standard_file)
        corrected_file = standard_file_path.parent / f'{standard_file_path.stem}_corrected.csv'

        # Rename columns by removing standard_space prefix
        column_mapping = {col: col.replace(f'{standard_space}_', '') 
                         for col in recon_standard.columns 
                         if f'{standard_space}_' in col}
        recon_standard = recon_standard.rename(columns=column_mapping)

        # Save to new file
        recon_standard.to_csv(corrected_file, index=False)

        return corrected_file

#%%
def run_pipeline(pre_implant_mri, 
                post_implant_ct, 
                ct_electrodes, 
                output_dir, 
                env_path=None, 
                freesurfer_dir=None,
                modules=['1', '2', '3', '4'], 
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
        freesurfer_dir (str/Path, optional): Path to FreeSurfer subjects directory
        modules (list): List of modules to run ['1', '2', '3', '4']
        skip_existing (bool): Skip processing if output files exist
        reg_type (str): Registration type ('gc', 'g', 'gc_noCTthereshold')
        qa_viewer (str): Quality assurance viewer type
    
    Returns:
        tuple: Paths to output files for modules 2, 3, and 4 (if run)
    """
    # Set project path
    project_path = Path(__file__).parent.parent
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
    file_locations_module2 = None
    file_locations_module3 = None
    file_locations_module4 = None
    
    if '1' in modules:
        print("Running Module 1...")
        recon.module1()
    
    if '2' in modules:
        print("Running Module 2...")
        file_locations_module2 = recon.module2(reg_type, skip_existing=skip_existing)
        
        print("Module 2 output files:")
        for name, path in file_locations_module2.items():
            print(f"{name}: {path}")
        
        recon.module2_QualityAssurance(file_locations_module2, qa_viewer)

    if '3' in modules:
        print("Running Module 3...")
        atlas = freesurfer_dir / 'mri' / 'aparc+aseg.mgz'
        atlas_lut = project_path / 'doc' / 'atlasLUT' / 'desikanKilliany.csv'
        file_locations_module3 = recon.module3(atlas, atlas_lut, diameter=2.5, skip_existing=skip_existing)
        
        print("Module 3 output file:")
        print(f"electrodes2ROI: {file_locations_module3}")
        
        recon.module3_QualityAssurance(file_locations_module3)

    if '4' in modules:
        print("Running Module 4...")
        file_locations_module4 = recon.module4(skip_existing=skip_existing)
        atlas_lut = project_path / 'doc' / 'atlasLUT' / 'desikanKilliany.csv'
        file_locations_module4_mni152 = recon.module4_snap_to_atlas(standard_space='mni152',
                                                          atlas='aparc+aseg.mgz',
                                                          atlas_lut=atlas_lut,
                                                          diameter=2.5)
        file_locations_module4_mni305 = recon.module4_snap_to_atlas(standard_space='mni305',
                                                          atlas='aparc+aseg.mgz',
                                                          atlas_lut=atlas_lut,
                                                          diameter=2.5)
        print(f"Module 4 output files:")
        print(f"MNI152: {file_locations_module4_mni152}")
        print(f"MNI305: {file_locations_module4_mni305}")

        # make a dictionary of the file locations for each module
        file_locations = {
            'module2': file_locations_module2,
            'module3': file_locations_module3,
            'module4': file_locations_module4,
            'module4_mni152': file_locations_module4_mni152,
            'module4_mni305': file_locations_module4_mni305
        }
    
    return file_locations

#%%
if __name__ == "__main__":
    # Example usage - replace these values with your actual file paths
    project_path = Path(__file__).parent.parent
   
    # Set paths for the selected subject
    pre_implant_mri = project_path / 'data' / 'ImagingData/NYU/sub-9740-10/ses-preimplant01/anat/sub-9740-10_ses-preimplant01_T1w.nii.gz'
    post_implant_ct = project_path / 'data' / 'ImagingData/NYU/sub-9740-10/ses-implant01/ct/sub-9740-10_ses-implant01_ct.nii.gz'
    ct_electrodes = project_path / 'data' / 'ImagingData/NYU/sub-9740-10/ses-implant01/ieeg/sub-9740-10_ses-implant01_t1_electrodes.tsv'
    output_dir = project_path / 'data' / 'output' / 'sub-9740-10'
    freesurfer_dir = None
   
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
        modules=['2'],
        skip_existing=True,
        reg_type='gc_noCTthereshold',  # Default registration type
        qa_viewer='niplot'  # Default viewer
    )
    
    print("Processing complete!")

# %%
