clc
clear
close all

%% Make ieeg_recon object

subject = ieeg_recon;

% by default the ieeg_recon object points to example data. This can be
% changed by setting the location of FILES: preImplantMRI, postImplantCT,
% postImplantCT_electrodes, and FOLDERS: output, and fslLoc in any object
% belonging to the class ieeg_recon

disp(subject)

subject.preImplantMRI = 'exampleData/sub-RID0922/ses-clinical01/anat/sub-RID0922_ses-clinical01_acq-3D_space-T00mri_T1w.nii.gz';
subject.postImplantCT = 'exampleData/sub-RID0922/ses-clinical01/ct/sub-RID0922_ses-clinical01_acq-3D_space-T01ct_ct.nii.gz';
subject.postImplantCT_electrodes = 'exampleData/sub-RID0922/ses-clinical01/ieeg/sub-RID0922_ses-clinical01_space-T01ct_desc-vox_electrodes.txt';
subject.output = 'exampleData/sub-RID0922/derivatives';
subject.fslLoc = '/usr/local/fsl/bin';

%% Run Module 1

subject.module1

%% Run Module 2

subject.module2