function [atlas,fs_name] = getAtlas(freesurfer_folder,atlasName)

FreeSurferColorLUT = readtable([freesurfer_folder '/FreeSurferColorLUT.txt']);
subject = 'fsaverage';

switch atlasName
    
    case 'subcortical'
        atlas = niftiread([freesurfer_folder '/' subject '/mri/aseg.nii.gz']);
        atlas = FreeSurferColorLUT(ismember(FreeSurferColorLUT.Var1,unique(atlas)),:);
        
        % Remove white matter, brain stem etc..
        atlas = removevars(atlas, {'Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13','Var14'});
        atlas([1:7,12:14,17,19,20,21:27,35,36,37:end],:) = [];
        
        atlas.Properties.VariableNames{1} = 'roiNum';
        atlas.Properties.VariableNames{2} = 'roi';
        atlas.Properties.VariableNames{3} = 'R';
        atlas.Properties.VariableNames{4} = 'G';
        atlas.Properties.VariableNames{5} = 'B';
        
        fs_name = 'aseg';
        
    case 'DKatlas'
        atlas = niftiread([freesurfer_folder '/' subject '/mri/aparc+aseg.nii.gz']);
        atlas = FreeSurferColorLUT(ismember(FreeSurferColorLUT.Var1,unique(atlas)),:);
        atlas([1:6,11:13,16,18,19,20:25,33,34,35:43,47,79,83],:) = [];
        atlas = removevars(atlas, {'Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13','Var14'});
        atlas.Properties.VariableNames{1} = 'roiNum';
        atlas.Properties.VariableNames{2} = 'roi';
        atlas.Properties.VariableNames{3} = 'R';
        atlas.Properties.VariableNames{4} = 'G';
        atlas.Properties.VariableNames{5} = 'B';

         fs_name = 'aparc+aseg';
        
    case 'Destrieux'
        atlas = niftiread([freesurfer_folder '/' subject '/mri/aparc.a2009s+aseg.nii.gz']);
        atlas = FreeSurferColorLUT(ismember(FreeSurferColorLUT.Var1,unique(atlas)),:);
        atlas([1:6,11:13,16,20:25,35:43,47,79,83],:) = [];
        atlas = removevars(atlas, {'Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13','Var14'});
        atlas.Properties.VariableNames{1} = 'roiNum';
        atlas.Properties.VariableNames{2} = 'roi';
        atlas.Properties.VariableNames{3} = 'R';
        atlas.Properties.VariableNames{4} = 'G';
        atlas.Properties.VariableNames{5} = 'B';

         fs_name = 'aparc.a2009s+aseg';
end


end
