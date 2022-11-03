function   cord_mm = getatlasROI(record_id,freesurfer_folder,atlas,atlas_fs_name)

subject = ['sub-RID' num2str(record_id,'%04.f')];
T1roi = atlas_fs_name;

T1.hdr = niftiinfo([freesurfer_folder '/' subject '/mri/' T1roi '.nii.gz']);
T1.data = niftiread([freesurfer_folder '/' subject '/mri/' T1roi '.nii.gz']);

cord_voxels = [];

for roi = 1:size(atlas,1)

    [CRS(:,1),CRS(:,2),CRS(:,3)]=ind2sub(size(T1.data),find(T1.data==atlas.roiNum(roi)));
    CRS(:,4) = repmat(atlas.roiNum(roi),size(CRS,1),1);
    cord_voxels = [cord_voxels;CRS];
    clear CRS
end

%% voxel coordinantes in mm space right anterior superior
cord_mm = T1.hdr.Transform.T' * [cord_voxels(:,1:3), ones(size(cord_voxels,1),1)]';
cord_mm = transpose(cord_mm);
cord_mm = [cord_mm(:,1:3), cord_voxels(:,4)];

end