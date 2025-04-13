#%%
import os
from pathlib import Path
from ieeg_recon.ieeg_recon import run_pipeline
import subprocess
import pandas as pd
from IPython import embed

#%%

def setup_subject_paths(project_path, subject_id):
    """
    Set up important file paths for a specific subject in the iEEG reconstruction pipeline.
    
    Parameters:
    -----------
    project_path : str or Path
        Base path to the project directory
    subject_id : str
        Subject ID (e.g., 'sub-RID0808')
        
    Returns:
    --------
    DataFrame
        DataFrame containing all important paths for the subject with each path type as a column
    """
    # Convert to Path object if string
    project_path = Path(project_path)
    subject_path = project_path / subject_id
    
    # Set up paths dictionary
    paths_dict = {
        # Pre-implant MRI from FreeSurfer
        't1_mgz': subject_path / 'derivatives' / 'freesurfer' / 'mri' / 'T1.mgz',
        't1': subject_path / 'derivatives' / 'freesurfer' / 'mri' / 'T1.nii.gz',
        
        # Output directory
        'output-dir': subject_path / 'derivatives', 
        
        # FreeSurfer directory
        'freesurfer-dir': subject_path / 'derivatives' / 'freesurfer'
    }
    
    # Find CT scan using glob (case-insensitive)
    ct_path = subject_path / 'ses-clinical01' / 'ct'
    ct_files = []
    for pattern in ['*CT.nii.gz', '*ct.nii.gz']:
        ct_files.extend(list(ct_path.rglob(pattern)))
    paths_dict['ct'] = ct_files[0] if ct_files else ''
    
    # Find electrode file using glob
    ieeg_path = subject_path / 'ses-clinical01' / 'ieeg'
    electrode_files = list(ieeg_path.rglob('*electrodes.txt'))
    paths_dict['elec'] = electrode_files[0] if electrode_files else ''
    
    # Convert T1.mgz to T1.nii.gz if needed
    if paths_dict['t1_mgz'].exists() and paths_dict['t1_mgz'].suffix == '.mgz' and not paths_dict['t1'].exists():
        paths_dict['t1'] = paths_dict['t1_mgz'].parent / 'T1.nii.gz'
        subprocess.run(['mri_convert', str(paths_dict['t1_mgz']), str(paths_dict['t1'])])
    
    # Create a DataFrame with each path type as a column and convert Path objects to strings
    paths_df = pd.DataFrame({
        'subject_id': [subject_id],
        't1_mgz': [str(paths_dict['t1_mgz'])],
        't1': [str(paths_dict['t1'])],
        'ct': [str(paths_dict['ct'])],
        'elec': [str(paths_dict['elec'])],
        'output-dir': [str(paths_dict['output-dir'])],
        'freesurfer-dir': [str(paths_dict['freesurfer-dir'])]
    })
    
    return paths_df

#%%

def get_all_subjects(project_path):
    """
    Find all subject directories in the project path.
    
    Parameters:
    -----------
    project_path : str or Path
        Base path to the project directory
        
    Returns:
    --------
    list
        List of subject IDs found in the project directory
    """
    project_path = Path(project_path)
    # Look for directories that match the pattern 'sub-RID*'
    subject_dirs = [d.name for d in project_path.glob('sub-RID*') if d.is_dir()]
    return subject_dirs

#%% 

def find_subject_data(project_path, subject_RID=None):
    """
    Find subject data for either a specific subject, multiple subjects, or all subjects in the project path.
    
    Parameters:
    -----------
    project_path : str or Path
        Base path to the project directory
    subject_RID : str, list, optional
        Subject ID (e.g., 'RID0808') or list of subject IDs. If provided, only those subjects' data is returned.
        If None, data for all subjects is returned.
        
    Returns:
    --------
    DataFrame
        DataFrame containing path information for the requested subject(s)
    """
    project_path = Path(project_path)
    df_all = pd.DataFrame()
    
    # Handle different input types for subject_RID
    if subject_RID is None:
        # Get all subjects if no specific subject is requested
        subjects = get_all_subjects(project_path)
    elif isinstance(subject_RID, list):
        # Process a list of subjects
        subjects = []
        for rid in subject_RID:
            # Check if each subject_RID already has the 'sub-' prefix
            if not str(rid).startswith('sub-'):
                subjects.append(f"sub-{rid}")
            else:
                subjects.append(rid)
    else:
        # Handle single subject as string
        if not subject_RID.startswith('sub-'):
            subject_id = f'sub-{subject_RID}'
        else:
            subject_id = subject_RID
        subjects = [subject_id]
    
    # Process each subject
    for subject_id in subjects:
        print(f"\nProcessing {subject_id}...")
        
        # Set up paths for this subject
        paths_df = setup_subject_paths(project_path, subject_id)
        
        # Check if required files exist
        required_files = ['t1', 'ct', 'elec']
        missing_files = [f for f in required_files if paths_df[f].values[0] == '']
        
        if missing_files:
            print(f"WARNING: Missing required files for {subject_id}:")
            for file in missing_files:
                print(f"  - {file}: {paths_df[file].values[0]}")
        else:
            print(f"All required files found for {subject_id}")

        df_all = pd.concat([df_all, paths_df])
    
    return df_all

#%%

if __name__ == "__main__":
    
    # Set the project path
    project_path = Path('/Users/nishant/Dropbox/Sinha/Lab/Research/epi_t3_iEEG/data/BIDS')
    
    # Example 1: Get data for all subjects
    all_subjects_df = find_subject_data(project_path)
    all_subjects_df.to_csv('all_subjects_paths.csv', index=False)
    
    # Example 2: Get data for a specific subject
    # specific_subject_df = find_subject_data(project_path, "RID0808")
    # specific_subject_df.to_csv('single_subject_paths.csv', index=False)
    
# %%
