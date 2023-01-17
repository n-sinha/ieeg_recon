classdef ieeg_recon

    %ieeg_recon is the electrode coregistration pipeline
    %   Detailed explanation goes here

    properties
        preImplantMRI
        postImplantCT
        postImplantCT_electrodes
        output
        fslLoc
        itksnap
    end

    methods

        function module1(obj)
            %module1: exports electrode coordinates of post implant CT in voxel and
            %native space. Outputs of this module goes in
            %output:ieeg_recon/module1 folder

            config_iEEGrecon;

            mkdir(obj.output, 'ieeg_recon/module1')

            % Export electrode cordinates in CT space in mm and vox
            elecCTvox = readtable(obj.postImplantCT_electrodes);
            writecell(elecCTvox{:, 1}, fullfile(obj.output, 'ieeg_recon/module1/electrode_names.txt'));
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

        end

        function module2(obj)
            %module2: Outputs of this module goes in
            %output:ieeg_recon/module2 folder

            config_iEEGrecon;
            mkdir(obj.output, 'ieeg_recon/module2');

            try

                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRImm.txt'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRIvox.txt'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRIvox.txt'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRI.nii.gz'));
                mustBeFile(fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRI_freesurferLUT.txt'));
               
            catch

                % remove negative values from CT image
                cmd = [obj.fslLoc '/fslmaths ' obj.postImplantCT ' -thr 0 ' ...
                    fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz')];

                system(cmd, "-echo");

                % FLIRT: Register CT to preimapnt T1 MRI
                cmd = [obj.fslLoc '/flirt' ...
                    ' -cost mutualinfo' ...
                    ' -dof 6' ...
                    ' -v' ...
                    ' -in ' fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -ref ' obj.preImplantMRI ...
                    ' -out ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -omat ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_flirt.txt')];
                %' -searchrx -180 180 -searchry -180 180  -searchrz -180 180'

                system(cmd, "-echo");

                % Greedy: Fine-tune flirt registered ct_to_mri to preimapnt T1 MRI
                cmd = [obj.itksnap '/greedy' ...
                    ' -d 3' ...
                    ' -i ' obj.preImplantMRI ...
                    ' ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -o ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_greedy.mat') ...
                    ' -a -dof 6' ...
                    ' -m NMI' ...
                    ' -ia-image-centers' ...
                    ' -n 100x100x0x0'];

                system(cmd, "-echo");

                % Convert greedy transform to FSL format
                cmd = [obj.itksnap '/c3d_affine_tool' ...
                    ' -ref ' obj.preImplantMRI ...
                    ' -src ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_greedy.mat') ...
                    ' -ras2fsl' ...
                    ' -o ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_greedy.txt')];

                system(cmd, "-echo");

                % Concatenate greedy and flirt transforms
                cmd = [obj.fslLoc '/convert_xfm' ...
                    ' -omat ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt') ...
                    ' -concat' ...
                    ' ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_greedy.txt') ...
                    ' ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_flirt.txt')];

                system(cmd, "-echo");

                % Apply combined registration to CT
                % replace flirt ct_to_mri with new ct_to_mri
                cmd = [obj.fslLoc '/flirt' ...
                    ' -in ' fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -ref ' obj.preImplantMRI ...
                    ' -init ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt') ...
                    ' -out ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -applyxfm'];

                system(cmd, "-echo");

                % Apply registration to electrode coordinates both in mm and vox
                cmd = [obj.fslLoc '/img2imgcoord' ...
                    ' -src ' fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -dest ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -xfm ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt') ...
                    ' -mm ' fullfile(obj.output, 'ieeg_recon/module1/electrodes_inCTmm.txt') ...
                    ' > ' fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRImm.txt')];

                system(cmd, "-echo")

                cmd = [obj.fslLoc '/img2imgcoord' ...
                    ' -src ' fullfile(obj.output, 'ieeg_recon/module2/ct_thresholded.nii.gz') ...
                    ' -dest ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri.nii.gz') ...
                    ' -xfm ' fullfile(obj.output, 'ieeg_recon/module2/ct_to_mri_xform.txt') ...
                    ' -vox ' fullfile(obj.output, 'ieeg_recon/module1/electrodes_inCTvox.txt') ...
                    ' > ' fullfile(obj.output, 'ieeg_recon/module2/electrodes_inMRIvox.txt')];

                system(cmd, "-echo");

                % Draw sphere of 2mm radius with electrode cordinate in
                % the center in the registered CT space

                % Create a blank nifti file in the space of registered ct_to_mri.nii.gz
                hdr = niftiinfo([obj.output '/ieeg_recon/module2/ct_to_mri.nii.gz']);
                data =  niftiread([obj.output '/ieeg_recon/module2/ct_to_mri.nii.gz']);
                data(data~=0) = 0;
                [dataVox(:,1), dataVox(:,2), dataVox(:,3)] = ind2sub(size(data), find(data==0));
                datamm = hdr.Transform.T' * [dataVox, ones(size(dataVox, 1), 1)]';
                datamm = transpose(datamm);
                datamm(:,4) = [];

                % Load electrode coordinates in mm
                electrodes_mm = importdata([obj.output '/ieeg_recon/module2/electrodes_inMRImm.txt']);
                electrodes_mm = electrodes_mm.data;

                % Export electrode labels for freesurfer
                electrode_freesurferLUT.index = transpose(1:size(electrodes_mm,1));
                electrode_freesurferLUT.names = importdata([obj.output '/ieeg_recon/module1/electrode_names.txt']);
                electrode_freesurferLUT.R = repmat(90, size(electrodes_mm,1),1);
                electrode_freesurferLUT.G = repmat(150, size(electrodes_mm,1),1);
                electrode_freesurferLUT.B = repmat(60, size(electrodes_mm,1),1);
                electrode_freesurferLUT.alpha = zeros(size(electrodes_mm,1),1);
                electrode_freesurferLUT = struct2table(electrode_freesurferLUT);
                writetable(electrode_freesurferLUT, [obj.output '/ieeg_recon/module2/electrodes_inMRI_freesurferLUT.txt'], ...
                    'Delimiter','space', 'WriteVariableNames', false);

                % Assign point for each electrode with a unique value
                [idx, dist_mm] = knnsearch(electrodes_mm,datamm);

                % only keep dist_mm <= 2 mm
                electrodesData = [find(dist_mm <= 2), idx(dist_mm <= 2), dist_mm(dist_mm <= 2)];
                [electrodesData(:,4), electrodesData(:,5), electrodesData(:,6)] = ind2sub(size(data), electrodesData(:,1));

                for n = 1:size(electrodesData,1)
                    data(electrodesData(n,4), electrodesData(n,5), electrodesData(n,6)) = electrodesData(n,2);
                end

                electrode_map = [obj.output '/ieeg_recon/module2/electrodes_inMRI.nii'];
                niftiwrite(data, electrode_map ,hdr, 'Compressed',true);
            end

            %% To DO: Export a sanpshot of CT over MRI to check registration quality

        end

        function electrodes2ROI = module3(obj, atlas, lookupTable)
            %module3: Outputs of this module goes in
            %output:ieeg_recon/module3 folder

            config_iEEGrecon;
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

            for r = 1:size(lut, 1)

                [vox(:, 1), vox(:, 2), vox(:, 3)] = ind2sub(size(data), find(data == lut.roiNum(r)));
                vox(:, 4) = repmat(lut.roiNum(r), size(vox, 1), 1);
                atlas_voxels = [atlas_voxels; vox];
                clear vox

            end

            %% Atlas roi in mm
            cord_mm = hdr.Transform.T' * [atlas_voxels(:, 1:3), ones(size(atlas_voxels, 1), 1)]';
            cord_mm = transpose(cord_mm);
            cord_mm = [cord_mm(:, 1:3), atlas_voxels(:, 4)];

            %% Map electrode contacts to all ROI

            [idx, dist_mm] = knnsearch(cord_mm(:, 1:3), electrodes_mm, 'K', 1);
            implant2roiNum = cord_mm(idx, 4);
            [~, idx] = ismember(implant2roiNum, lut.roiNum);
            implant2roi = lut.roi(idx);

            % if an electrode is not within 2.5mm of any ROI it is in the white matter
            % or outside the brain
            implant2roiNum(dist_mm >= 2.6) = nan;
            [implant2roi{dist_mm >= 2.6, :}] = deal('');

            electrodes2ROI = table(labels, ...
                electrodes_mm(:, 1), electrodes_mm(:, 2), electrodes_mm(:, 3), ...
                electrodes_vox(:, 1), electrodes_vox(:, 2), electrodes_vox(:, 3), ...
                implant2roi, implant2roiNum);

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
