#!/bin/bash

# Load required modules
module load ANTs/2.3.5
module load fsl/6.0.3
module load ITK/5.2.1 vtk/9.2.6 greedy
module load openblas/0.3.27
module load freesurfer/7.1.0
source /appl/freesurfer-7.1.0/SetUpFreeSurfer.sh
export FS_LICENSE=/project/davis_group_1/nishants/license.txt
export SURFER_FRONTDOOR=1

rid=$1

cd /project/davis_group_1/nishants/ieeg_recon
source .venv/bin/activate

# Run the Python script with input arguments
python run_ieeg_recon.py \
    --t1 /project/davis_group_1/nishants/ieeg_recon/data/BIDS/${rid}/derivatives/freesurfer/mri/T1.nii.gz \
    --ct /project/davis_group_1/nishants/ieeg_recon/data/BIDS/${rid}/derivatives/freesurfer/mri/T1.nii.gz \
    --elec /project/davis_group_1/nishants/ieeg_recon/data/BIDS/${rid}/derivatives/freesurfer/mri/T1.nii.gz \
    --output-dir /project/davis_group_1/nishants/ieeg_recon/data/BIDS/${rid}/derivatives \
    --freesurfer-dir /project/davis_group_1/nishants/ieeg_recon/data/BIDS/${rid}/derivatives/freesurfer \
    --modules 4 \