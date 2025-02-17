#%% 
import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
from pathlib import Path
from scipy.spatial import cKDTree
import argparse
from dotenv import load_dotenv
from nibabel.freesurfer.io import read_geometry  # For reading FreeSurfer surfaces
import trimesh  # For point-in-polyhedron tests

#%% 
class IEEGRecon:
    def __init__(self, pre_implant_mri, post_implant_ct, ct_electrodes, output_dir, env_path=None):
        """
        Initialize IEEGRecon with required paths
        
        Args:
            pre_implant_mri (str): Path to pre-implant MRI
            post_implant_ct (str): Path to post-implant CT
            ct_electrodes (str): Path to electrode coordinates CSV
            output_dir (str): Output directory path (required)
            env_path (str, optional): Path to .env file
        """
        # Set main properties
        self.preImplantMRI = pre_implant_mri
        self.postImplantCT = post_implant_ct
        self.postImplantCT_electrodes = ct_electrodes
        self.output = Path(output_dir)  # Convert to Path object
        self.root_dir = Path(__file__).parent.parent
        
        # Setup environment variables and paths
        self._setup_environment(env_path)

    def _setup_environment(self, env_path=None):
        """Setup environment variables and paths from .env file"""
        if env_path is None:
            env_path = Path(__file__).parent / '.env'
        
        try:
            # Load environment variables from .env file
            load_dotenv(env_path)
            
            # Set paths from environment variables
            self.fslLoc = os.getenv('FSL_DIR')
            self.itksnap = os.getenv('ITKSNAP_DIR')
            self.freeSurfer = os.getenv('FREESURFER_HOME')
            self.freeSurferDir = os.getenv('SUBJECTS_DIR')
            
            # Ensure required environment variables are set
            if not all([self.fslLoc, self.itksnap, self.freeSurfer, self.freeSurferDir]):
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
            imageviewer (str): Type of viewer to use - 'freeview_snapshot', 'freeview', or 'itksnap'
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
            else:
                raise ValueError(f"Unknown imageviewer option: {imageviewer}")

    def module3(self, atlas, lookup_table, diameter=2.5, skip_existing=False):
        """
        Module3: Map electrodes to brain regions using provided atlas
        
        Args:
            atlas (str/Path): Path to atlas NIFTI file
            lookup_table (str/Path): Path to lookup table CSV/txt file
            diameter (float): Maximum distance in mm for electrode-to-ROI mapping (default: 2.5)
            skip_existing (bool): If True, skip processing if output files exist
        
        Returns:
            str: Path to output electrodes2ROI CSV file
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon/module3'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define output file location
        output_file = output_dir / 'electrodes2ROI.csv'
        
        # Check if file exists and skip if requested
        if skip_existing and output_file.exists():
            print("Module 3 output already exists, skipping...")
            return str(output_file)
        
        # Verify input files exist
        for file in [atlas, lookup_table]:
            if not Path(file).exists():
                raise FileNotFoundError(f"Input file not found: {file}")

        # Load electrode coordinates and names
        electrodes_mm = np.loadtxt(
            Path(self.output) / 'ieeg_recon/module2/electrodes_inMRImm.txt',
            skiprows=1  # Skip header if present
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
        lut = pd.read_csv(lookup_table, sep=None, engine='python')  # Automatically detect separator

        # Get atlas ROI coordinates in voxels
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
        implant2roiNum[dist_mm > diameter] = np.nan
        
        # Load FreeSurfer surfaces
        lh_pial, _ = read_geometry(Path(self.freeSurferDir) / 'surf/lh.pial')
        lh_white, _ = read_geometry(Path(self.freeSurferDir) / 'surf/lh.white')
        rh_pial, _ = read_geometry(Path(self.freeSurferDir) / 'surf/rh.pial')
        rh_white, _ = read_geometry(Path(self.freeSurferDir) / 'surf/rh.white')
        
        # Convert surfaces to trimesh format
        lh_pial_mesh = trimesh.Trimesh(vertices=lh_pial, faces=_)
        rh_pial_mesh = trimesh.Trimesh(vertices=rh_pial, faces=_)
        lh_white_mesh = trimesh.Trimesh(vertices=lh_white, faces=_)
        rh_white_mesh = trimesh.Trimesh(vertices=rh_white, faces=_)
        
        # Find contacts outside brain (not in either pial surface)
        outside_mask = ~(
            lh_pial_mesh.contains(electrodes_mm[dist_mm > diameter]) | 
            rh_pial_mesh.contains(electrodes_mm[dist_mm > diameter])
        )
        outside_indices = np.where(dist_mm > diameter)[0][outside_mask]
        implant2roi.iloc[outside_indices] = 'outside-brain'
        
        # Find white matter contacts (inside white surface)
        wm_indices = np.where(dist_mm > diameter)[0][~outside_mask]
        implant2roi.iloc[wm_indices] = 'white-matter'
        
        # Verify white matter contacts
        wm_contacts = electrodes_mm[wm_indices]
        in_white = (
            lh_white_mesh.contains(wm_contacts) | 
            rh_white_mesh.contains(wm_contacts)
        )
        
        if np.all(in_white):
            print('white-matter contacts correctly assigned')
        else:
            print(f"check white-matter contacts in: {self.output}")
        
        # Create output dataframe
        electrodes2ROI = pd.DataFrame({
            'labels': labels,
            'mm_x': electrodes_mm[:, 0],
            'mm_y': electrodes_mm[:, 1],
            'mm_z': electrodes_mm[:, 2],
            'vox_x': electrodes_vox[:, 0],
            'vox_y': electrodes_vox[:, 1],
            'vox_z': electrodes_vox[:, 2],
            'roi': implant2roi,
            'roiNum': implant2roiNum
        })
        
        # Save output
        electrodes2ROI.to_csv(output_file, index=False)
        
        return str(output_file)

def run_pipeline(pre_implant_mri, post_implant_ct, ct_electrodes, output_dir, env_path=None,
                modules=['1', '2'], skip_existing=False, reg_type='gc_noCTthereshold', qa_viewer='freeview'):
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
        env_path=env_path
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
        atlas = '/Users/nishant/Dropbox/Sinha/Lab/Research/iEEG_recon_local/data/sub-RID0031/derivatives/freesurfer/mri/aparc+aseg.nii.gz'
        lookup_table = project_path / 'doc' / 'atlasLUT' / 'desikanKilliany.csv'
        recon.module3(atlas, lookup_table, diameter=2.5)
    
    return file_locations

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
        modules=['3'],  # Run both modules by default
        skip_existing=True,  # Don't skip existing files by default
        reg_type='gc_noCTthereshold',  # Default registration type
        qa_viewer='freeview'  # Default viewer
    )
    
    print("Processing complete!")

# %%
