clc
clear
close all

%% Make ieeg_recon object

subject = ieeg_recon;

% by default the ieeg_recon object points to example data. This can be
% changed by setting the location of FILES: preImplantMRI, postImplantCT,
% postImplantCT_electrodes, and FOLDERS: output, and fslLoc in any object
% belonging to the class ieeg_recon

subject.preImplantMRI = '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/BIDS/derivatives/t3_freesurfer/sub-RID0031/mri/T1.nii.gz';
subject.postImplantCT = '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/BIDS/derivatives/t3_freesurfer/sub-RID0031/CT/native/CT.nii.gz';
subject.postImplantCT_electrodes = '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/BIDS/derivatives/t3_freesurfer/sub-RID0031/CT/native/electrodes_inCTcrs.txt';
subject.output = '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/BIDS/derivatives/t3_freesurfer/sub-RID0031';
subject.fslLoc = '/usr/local/fsl/bin';

%% Run Module 1

subject.module1

%% Run Module 2

subject.module2

%% Run Module 3

atlas = '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/BIDS/derivatives/t3_freesurfer/sub-RID0031/mri/aparc+aseg.nii.gz';
lookupTable = 'atlas_lookuptable/desikanKilliany.csv';
 
subject.module3(atlas, lookupTable)
