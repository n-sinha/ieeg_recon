clc
clear
close all
warning off
addpath(genpath('external_lib'));

%% Load atlas and keep only the ROI of interest

freesurfer_folder = 'Data/freesurfer';
atlasName = 'DKatlas'; % select from 'subcortical', 'DKatlas', or 'Destrieux' atlas

[atlas,atlas_fs_name] = getAtlas(freesurfer_folder,atlasName);

%region_of_interest_lh = atlas.roi(15:48); % for all ROIs in left side
region_of_interest_lh = {'ctx-lh-insula'};

%region_of_interest_rh = atlas.roi(49:82);
region_of_interest_rh = {'ctx-rh-insula'}; % for all ROIs on right side

%% Make plots for every subject

record_id = [31;32;37;50;51;59;89;101;102;117;139;...
    143;179;190;213;278;309;320;365;420;...
    440;454;459;490;502;508;520;522;529;...
    536;566;572;583;595;596;646;648;652;...
    656;658;679;785];


for sub = 1%:size(record_id,1)

   
    subject = ['sub-RID' num2str(record_id(sub),'%04.f')];
    load(['Data/iEEG_coord/' subject '_cord.mat']);

    [lpv, lpf] = read_surf([freesurfer_folder '/' subject '/surf/lh.pial']);
    [rpv, rpf] = read_surf([freesurfer_folder '/' subject '/surf/rh.pial']);

    %% Figure 1:  Plot pial surface of subject and all electrodes grouped by region
    figure;

    hold on
    hl = trisurf(lpf+1,lpv(:,1),lpv(:,2),lpv(:,3));
    hl.LineStyle = 'none';
    hl.FaceAlpha = 0.1;
    hl.FaceColor = [0.7 0.7 0.7];
    hr = trisurf(rpf+1,rpv(:,1),rpv(:,2),rpv(:,3));
    hr.LineStyle = 'none';
    hr.FaceAlpha = 0.1;
    hr.FaceColor = [0.9 0.7 0.8];

    groups = unique(cord.roiNum);
    for i = 1:length(groups)

        idx = cord.roiNum==groups(i);
        roiName(i,:) = unique(cord.roiName(idx));
        cordTemp = cord(idx,:);
        scatter3(cordTemp.cord_surf(:,1),cordTemp.cord_surf(:,2),...
            cordTemp.cord_surf(:,3),100*ones(size(cordTemp,1),1),cordTemp.roiColor,'filled')

    end
    legend(['lh'; 'rh'; roiName],'Location','westoutside','FontSize',16)
    hold off
    title(subject)

    %% Figure 2:  Plot pial surface of subject and all electrodes in a specific region

    figure;

    hold on
    hl = trisurf(lpf+1,lpv(:,1),lpv(:,2),lpv(:,3));
    hl.LineStyle = 'none';
    hl.FaceAlpha = 0.1;
    hl.FaceColor = [0.7 0.7 0.7];
    hr = trisurf(rpf+1,rpv(:,1),rpv(:,2),rpv(:,3));
    hr.LineStyle = 'none';
    hr.FaceAlpha = 0.1;
    hr.FaceColor = [0.7 0.7 0.7];

    % plot left side roi
    atlas_roi_lh = atlas(contains(atlas.roi,region_of_interest_lh),:);
    for roi = 1:size(atlas_roi_lh,1)
        [lhROIv, lhROIf] = getROISurfaceCord(freesurfer_folder,subject,atlas_roi_lh(roi,:),lpv, lpf,'lh');

        hl_roi = trisurf(lhROIf+1,lhROIv(:,1),lhROIv(:,2),lhROIv(:,3));
        hl_roi.LineStyle = 'none';
        %hl_roi.FaceAlpha = 0.1;
        hl_roi.FaceColor = atlas_roi_lh{roi,3:5}/255;
    end

    % plot right side roi
    atlas_roi_rh = atlas(contains(atlas.roi,region_of_interest_rh),:);
    for roi = 1:size(atlas_roi_rh,1)
        [rhROIv, rhROIf] = getROISurfaceCord(freesurfer_folder,subject,atlas_roi_rh(roi,:),rpv, rpf,'rh');

        hr_roi = trisurf(rhROIf+1,rhROIv(:,1),rhROIv(:,2),rhROIv(:,3));
        hr_roi.LineStyle = 'none';
        %hr_roi.FaceAlpha = 0.1;
        hr_roi.FaceColor = atlas_roi_rh{roi,3:5}/255;
    end

   scatter3(cord.cord_surf(:,1),cord.cord_surf(:,2),...
        cord.cord_surf(:,3),100*ones(size(cord,1),1),cord.roiColor,'filled','MarkerEdgeColor','k');
    hold off
    title(subject)

end