%% Set external library in path

addpath(genpath('external_lib'))

%% FSL Setup

setenv( 'FSLDIR', '/usr/local/fsl' );
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

%% ITK Snap Setup

setenv('ITKSNAPDIR', '/Applications/ITK-SNAP.app/Contents/bin');
itksnapdir = getenv('ITKSNAPDIR');
itksnapmpath = sprintf('%s',itksnapdir);
path(path,itksnapmpath)
clear itksnapdir itksnapmpath;

%% Freesurfer setup

setenv( 'FREESURFER_HOME', '/Applications/freesurfer/7.3.2/' );
setenv('SUBJECTS_DIR', '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/');
FREESURFER_HOME = getenv('FREESURFER_HOME');
freesurferdirmpath = sprintf('%s/SetUpFreeSurfer.sh',FREESURFER_HOME);
system(['sh ' freesurferdirmpath],'-echo');
clear freesurferdirmpath FREESURFER_HOME;
