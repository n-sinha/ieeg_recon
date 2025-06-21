#!/usr/bin/env python3
"""
Script to copy iEEG reconstruction data and related files from source to target BIDS directory.

This script:
1. Finds all subjects in the source BIDS directory that have 'derivatives/ieeg_recon' directory
2. Copies the entire 'ieeg_recon' directory structure to the target location
3. Copies specific FreeSurfer MRI files (T1.nii.gz, T1.mgz, brain.mgz) for the same subjects

Author: AI Assistant
Date: 2024
"""

import os
import shutil
import glob
from pathlib import Path
import logging

# Set up logging to track the copying process
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('copy_ieeg_recon.log'),
        logging.StreamHandler()
    ]
)

def find_subjects_with_ieeg_recon(source_bids_dir):
    """
    Find all subjects that have the 'derivatives/ieeg_recon' directory.
    
    Args:
        source_bids_dir (str): Path to the source BIDS directory
        
    Returns:
        list: List of subject IDs that have ieeg_recon directory
    """
    subjects_with_ieeg_recon = []
    
    # Get all subject directories (they start with 'sub-')
    subject_dirs = glob.glob(os.path.join(source_bids_dir, "sub-*"))
    
    for subject_dir in subject_dirs:
        subject_id = os.path.basename(subject_dir)
        ieeg_recon_path = os.path.join(subject_dir, "derivatives", "ieeg_recon")
        
        # Check if the ieeg_recon directory exists
        if os.path.exists(ieeg_recon_path) and os.path.isdir(ieeg_recon_path):
            subjects_with_ieeg_recon.append(subject_id)
            logging.info(f"Found ieeg_recon directory for {subject_id}")
    
    return subjects_with_ieeg_recon

def copy_ieeg_recon_directory(source_subject_dir, target_subject_dir):
    """
    Copy the entire ieeg_recon directory structure from source to target.
    
    Args:
        source_subject_dir (str): Source subject directory path
        target_subject_dir (str): Target subject directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    source_ieeg_recon = os.path.join(source_subject_dir, "derivatives", "ieeg_recon")
    target_ieeg_recon = os.path.join(target_subject_dir, "derivatives", "ieeg_recon")
    
    try:
        # Create the target derivatives directory if it doesn't exist
        os.makedirs(os.path.dirname(target_ieeg_recon), exist_ok=True)
        
        # Copy the entire ieeg_recon directory
        if os.path.exists(target_ieeg_recon):
            logging.warning(f"Target directory {target_ieeg_recon} already exists. Skipping copy.")
            return True
        
        shutil.copytree(source_ieeg_recon, target_ieeg_recon)
        logging.info(f"Successfully copied ieeg_recon directory for {os.path.basename(source_subject_dir)}")
        return True
        
    except Exception as e:
        logging.error(f"Error copying ieeg_recon directory for {os.path.basename(source_subject_dir)}: {str(e)}")
        return False

def copy_freesurfer_files(source_subject_dir, target_subject_dir):
    """
    Copy specific FreeSurfer MRI files from source to target.
    
    Args:
        source_subject_dir (str): Source subject directory path
        target_subject_dir (str): Target subject directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    # Define the files to copy
    files_to_copy = [
        "derivatives/freesurfer/mri/T1.nii.gz",
        "derivatives/freesurfer/mri/T1.mgz", 
        "derivatives/freesurfer/mri/brain.mgz"
    ]
    
    success = True
    
    for file_rel_path in files_to_copy:
        source_file = os.path.join(source_subject_dir, file_rel_path)
        target_file = os.path.join(target_subject_dir, file_rel_path)
        
        # Check if source file exists
        if not os.path.exists(source_file):
            logging.warning(f"Source file not found: {source_file}")
            success = False
            continue
        
        try:
            # Create the target directory if it doesn't exist
            os.makedirs(os.path.dirname(target_file), exist_ok=True)
            
            # Copy the file
            shutil.copy2(source_file, target_file)
            logging.info(f"Successfully copied: {file_rel_path}")
            
        except Exception as e:
            logging.error(f"Error copying {file_rel_path}: {str(e)}")
            success = False
    
    return success

def main():
    """
    Main function to orchestrate the copying process.
    """
    # Define source and target directories
    source_bids_dir = "/Users/nishant/Dropbox/Sinha/Lab/Research/epi_t3_iEEG/data/BIDS"
    target_bids_dir = "/Users/nishant/Dropbox/Sinha/Lab/Research/iEEG_recon_local/data/BIDS"
    
    logging.info("Starting iEEG reconstruction data copy process")
    logging.info(f"Source BIDS directory: {source_bids_dir}")
    logging.info(f"Target BIDS directory: {target_bids_dir}")
    
    # Check if source directory exists
    if not os.path.exists(source_bids_dir):
        logging.error(f"Source BIDS directory does not exist: {source_bids_dir}")
        return
    
    # Create target BIDS directory if it doesn't exist
    os.makedirs(target_bids_dir, exist_ok=True)
    logging.info(f"Created target BIDS directory: {target_bids_dir}")
    
    # Find subjects with ieeg_recon directory
    subjects_with_ieeg_recon = find_subjects_with_ieeg_recon(source_bids_dir)
    
    if not subjects_with_ieeg_recon:
        logging.warning("No subjects found with ieeg_recon directory")
        return
    
    logging.info(f"Found {len(subjects_with_ieeg_recon)} subjects with ieeg_recon directory")
    
    # Process each subject
    successful_copies = 0
    failed_copies = 0
    
    for subject_id in subjects_with_ieeg_recon:
        logging.info(f"Processing subject: {subject_id}")
        
        source_subject_dir = os.path.join(source_bids_dir, subject_id)
        target_subject_dir = os.path.join(target_bids_dir, subject_id)
        
        # Create target subject directory
        os.makedirs(target_subject_dir, exist_ok=True)
        
        # Copy ieeg_recon directory
        ieeg_recon_success = copy_ieeg_recon_directory(source_subject_dir, target_subject_dir)
        
        # Copy FreeSurfer files
        freesurfer_success = copy_freesurfer_files(source_subject_dir, target_subject_dir)
        
        if ieeg_recon_success and freesurfer_success:
            successful_copies += 1
            logging.info(f"Successfully processed all files for {subject_id}")
        else:
            failed_copies += 1
            logging.error(f"Failed to process some files for {subject_id}")
    
    # Summary
    logging.info("=" * 50)
    logging.info("COPY PROCESS SUMMARY")
    logging.info("=" * 50)
    logging.info(f"Total subjects with ieeg_recon: {len(subjects_with_ieeg_recon)}")
    logging.info(f"Successfully processed: {successful_copies}")
    logging.info(f"Failed to process: {failed_copies}")
    logging.info("=" * 50)
    
    if successful_copies > 0:
        logging.info("Copy process completed successfully!")
    else:
        logging.error("No subjects were successfully processed!")

if __name__ == "__main__":
    main() 