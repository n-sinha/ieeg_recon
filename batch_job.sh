#!/bin/bash

# Path to the CSV file
csv_file="all_subjects_paths_process2.csv"

# Skip the header line and process each subject
tail -n +2 "$csv_file" | while IFS=, read -r subject_id t1_mgz t1 ct elec output_dir freesurfer_dir; do
    # Print which subject is being processed
    echo "=========================================="
    echo "Processing subject: $subject_id"
    echo "=========================================="
    
    # Run the ieeg_recon.py script with parameters from the CSV
    # Using python to execute the script
    python /Users/nishant/Dropbox/Sinha/Lab/Research/iEEG_recon_local/run_ieeg_recon.py \
        --t1 "$t1" \
        --ct "$ct" \
        --elec "$elec" \
        --output-dir "$output_dir" \
        --freesurfer-dir "$freesurfer_dir" \
        --qa-viewer 'niplot' \
        --reg-type 'gc_noCTthereshold' \
        --modules '4'
    
    # Print completion message
    echo "Completed processing $subject_id"
    echo "==========================================\n"
done

echo "All subjects have been processed!"
