classdef ieeg_recon

    %ieeg_recon is the electrode coregistration pipeline
    %   Detailed explanation goes here

    properties
        preImplantMRI {mustBeFile} = 'exampleData/sub-RID0922/ses-clinical01/anat/sub-RID0922_ses-clinical01_acq-3D_space-T00mri_T1w.nii.gz'
        postImplantCT {mustBeFile} = 'exampleData/sub-RID0922/ses-clinical01/ct/sub-RID0922_ses-clinical01_acq-3D_space-T01ct_ct.nii.gz'
        postImplantCT_electrodes {mustBeFile} = 'exampleData/sub-RID0922/ses-clinical01/ieeg/sub-RID0922_ses-clinical01_space-T01ct_desc-vox_electrodes.txt'
        output {mustBeFolder} = 'exampleData/sub-RID0922/derivatives'
        fslLoc {mustBeFolder} = '/usr/local/fsl/bin';
    end

    methods

        function module1(obj)
            %module1: exports electrode coordinates of post implant CT in voxel and
            %native space. Outputs of this module goes in
            %output:ieeg_recon/module1 folder

            config;

            mkdir(obj.output, 'ieeg_recon/module1')

            % Export electrode cordinates in CT space in mm and vox
            elecCTvox = readtable(obj.postImplantCT_electrodes);
            writecell(elecCTvox{:,1}, fullfile(obj.output, 'ieeg_recon/module1/electrode_names.txt'));
            writematrix(elecCTvox{:, 2:4}, fullfile(obj.output, 'ieeg_recon/module1/electrodes_inCTvox.txt'), 'Delimiter', 'space');

            CT.hdr = niftiinfo(obj.postImplantCT);
            CT.data = niftiread(obj.postImplantCT);

            threshold_to_removeSkull = 99.95; % percentile
            [CT.vox(:, 1), CT.vox(:, 2), CT.vox(:, 3)] = ind2sub(size(CT.data), ...
                find(CT.data > prctile(CT.data(:), threshold_to_removeSkull)));

            CT.mm = CT.hdr.Transform.T' * [CT.vox, ones(size(CT.vox, 1), 1)]';
            CT.mm = transpose(CT.mm);

            elecCTmm = transpose(CT.hdr.Transform.T' * [elecCTvox{:, 2:4}, ones(size(elecCTvox, 1), 1)]');
            writematrix(elecCTmm(:, 1:3), fullfile(obj.output, 'ieeg_recon/module1/electrodes_inCTmm.txt'), 'Delimiter', 'space')

            % export figure: electrode coordinates in mm ploted on CT native
            figure;
            hold on
            subplot(2, 3, 1);
            h1 = scatter3(CT.mm(:, 1), CT.mm(:, 2), CT.mm(:, 3), ...
                'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
            hold on
            h2 = scatter3(elecCTmm(:, 1), elecCTmm(:, 2), elecCTmm(:, 3), 'r', 'filled');
            hold off
            plotObjs = [h1, h2];
            plotinMultipleViews(plotObjs)
            sgtitle('Electrodes in CT mm space');
            fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [16 10] * 1.5);
            print(gcf, '-dpdf', '-r100', ...
                fullfile(obj.output, 'ieeg_recon/module1/check_quality.pdf'));
        end

        function module2(obj)
            %module2: Outputs of this module goes in
            %output:ieeg_recon/module2 folder

            config;

            mkdir(obj.output, 'ieeg_recon/module2');

            %% remove negative values from CT image
            try

                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz'));

            catch
                cmd = [obj.fslLoc '/fslmaths ' obj.postImplantCT ' -thr 0 ' ...
                    fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz')];

                system(cmd, "-echo");
            end

            %% Register CT to preimapnt T1 MRI
            try

                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt'));

            catch

                cmd = [obj.fslLoc '/flirt' ...
                    ' -cost mutualinfo' ...
                    ' -dof 6' ...
                    ' -v' ...
                    ' -in ' fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -ref ' obj.preImplantMRI ...
                    ' -searchrx -180 180' ...
                    ' -searchry -180 180' ...
                    ' -searchrz -180 180' ...
                    ' -out ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -omat ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt')];

                system(cmd, "-echo");

            end

            %% Apply registration to electrode coordinates both in mm and vox

            try

                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRImm.txt'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRIvox.txt'));

            catch

                cmd = [obj.fslLoc '/img2imgcoord' ...
                    ' -src '  fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -dest ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -xfm ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt') ...
                    ' -mm ' fullfile(obj.output, 'ieeg_recon/module1/electrodes_inCTmm.txt') ...
                    ' > ' fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRImm.txt')];

                system(cmd, "-echo")

                cmd = [obj.fslLoc '/img2imgcoord' ...
                    ' -src '  fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -dest ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -xfm ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt') ...
                    ' -vox ' fullfile(obj.output, 'ieeg_recon/module1/electrodes_inCTvox.txt') ...
                    ' > ' fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRIvox.txt')];

                system(cmd, "-echo");

            end

            %% Write electrode cordinates to nifti file

            % To DO

            %% Export a sanpshot of CT over MRI to check registration quality

            % To Do

        end

        function electrodes2ROI = module3(obj, atlas, lookupTable)
            %module3: Outputs of this module goes in
            %output:ieeg_recon/module3 folder
             
            config;
            mkdir(obj.output, 'ieeg_recon/module3');
            
            mustBeFile(atlas);
            mustBeFile(lookupTable);

            %% Load electrode coordinates in mm, and read atlas in native space from nifti file

            electrodes_mm = importdata([obj.output '/ieeg_recon/module2/electrodes_inMRImm.txt']);
            electrodes_mm = electrodes_mm.data;

            electrodes_vox = importdata([obj.output '/ieeg_recon/module2/electrodes_inMRIvox.txt']);
            electrodes_vox = electrodes_vox.data;

            labels = importdata([obj.output '/ieeg_recon/module1/electrode_names.txt']);

            hdr = niftiinfo(atlas);
            data = niftiread(atlas);
            lut = readtable(lookupTable);

            %% Atlas roi in voxels
            atlas_voxels = [];

            for r = 1:size(lut,1)

                [vox(:,1),vox(:,2),vox(:,3)]=ind2sub(size(data),find(data==lut.roiNum(r)));
                vox(:,4) = repmat(lut.roiNum(r),size(vox,1),1);
                atlas_voxels = [atlas_voxels;vox];
                clear vox

            end

            %% Atlas roi in mm
            cord_mm = hdr.Transform.T' * [atlas_voxels(:,1:3), ones(size(atlas_voxels,1),1)]';
            cord_mm = transpose(cord_mm);
            cord_mm = [cord_mm(:,1:3), atlas_voxels(:,4)];

            %% Map electrode contacts to all ROI

            [idx,dist_mm] = knnsearch(cord_mm(:,1:3),electrodes_mm,'K',1);
            implant2roiNum = cord_mm(idx,4);
            [~,idx] = ismember(implant2roiNum,lut.roiNum);
            implant2roi = lut.roi(idx);

            % if an electrode is not within 2.5mm of any ROI it is in the white matter
            % or outside the brain
            implant2roiNum(dist_mm>=2.6) = nan;
            [implant2roi{dist_mm>=2.6,:}] = deal('');


            electrodes2ROI = table(labels,...
                electrodes_mm(:,1),electrodes_mm(:,2),electrodes_mm(:,3),...
                electrodes_vox(:,1), electrodes_vox(:,2), electrodes_vox(:,3), ...
                implant2roi,implant2roiNum);

            electrodes2ROI.Properties.VariableNames(2) = "mm_x";
            electrodes2ROI.Properties.VariableNames(3) = "mm_y";
            electrodes2ROI.Properties.VariableNames(4) = "mm_z";
            electrodes2ROI.Properties.VariableNames(5) = "vox_x";
            electrodes2ROI.Properties.VariableNames(6) = "vox_y";
            electrodes2ROI.Properties.VariableNames(7) = "vox_z";
            electrodes2ROI.Properties.VariableNames(8) = "roi";
            electrodes2ROI.Properties.VariableNames(9) = "roiNum";

            writetable(electrodes2ROI, [obj.output, '/ieeg_recon/module3/electrodes2ROI.csv'])

        end

    end

end