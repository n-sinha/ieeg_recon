clc
clear
close all
warning off
addpath(genpath('external_lib'));

%% Load atlas and keep only the ROI of interest

freesurfer_folder = '/Users/nishant/Dropbox/Database/UPenn/Epilepsy/BIDS/derivatives/t3_freesurfer';
atlasName = 'DKatlas'; % select from 'subcortical', 'DKatlas', or 'Destrieux' atlas

[atlas,atlas_fs_name] = getAtlas(freesurfer_folder,atlasName);
region_of_interest = {'ctx-lh-insula','ctx-rh-insula'};

%% Load metadata exported from redcap

record_id = [31;32;37;50;51;59;89;101;102;117;139;...
    143;179;190;213;278;309;320;365;420;...
    440;454;459;490;502;508;520;522;529;...
    536;566;572;583;595;596;646;648;652;...
    656;658;679;785];

for sub = 1%:size(record_id,1)

    subject = ['sub-RID' num2str(record_id(sub),'%04.f')];
    load(['Data/iEEG_coord/' subject '_cord.mat']);

    cordROI = getatlasROI(record_id(sub),freesurfer_folder,atlas,atlas_fs_name);

    %% Map electrodes to all ROI

    [idx,dist_mm] = knnsearch(cordROI(:,1:3),cord.cord_mm,'K',1);
    implant2roi = cordROI(idx,4);
    % if an electrode is not within 2.5mm of any ROI it is in the white matter
    % or outside the brain
    implant2roi(dist_mm>2.5) = nan;

    % Check coverage in region_of_interest
    atlas_roi = atlas(contains(atlas.roi,region_of_interest),:);

    figure;
    hold on;
    scatter3(cord.cord_mm(:,1),cord.cord_mm(:,2),...
        cord.cord_mm(:,3),100*ones(size(cord,1),1),'k','MarkerEdgeColor','k');

    for roi = 1:size(atlas_roi)

        idx_ieeg = ismember(implant2roi,atlas_roi.roiNum(roi));
        atlas_roi.iEEGcoverage(roi) = sum(idx_ieeg);

        scatter3(cord.cord_mm(idx_ieeg,1),cord.cord_mm(idx_ieeg,2),...
            cord.cord_mm(idx_ieeg,3),100*ones(sum(idx_ieeg),1),...
            'r','filled','MarkerEdgeColor','k');

        idx_roi = ismember(cordROI(:,4),atlas_roi.roiNum(roi));
        scatter3(cordROI(idx_roi,1),cordROI(idx_roi,2),cordROI(idx_roi,3));

    end
   title(subject)
   hold off


end
