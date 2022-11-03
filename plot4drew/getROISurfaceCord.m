function [roiVert, roiFace] = getROISurfaceCord(freesurfer_folder,subject,atlas_roi,allvertices,allfaces,side)

[annot.vertices, annot.label, annot.colortable] = ...
    read_annotation([freesurfer_folder '/' subject '/label/' side '.aparc.annot']);


temp1 = and(and(ismember(annot.colortable.table(:,1), atlas_roi.R),...
    ismember(annot.colortable.table(:,2), atlas_roi.G)),...
    ismember(annot.colortable.table(:,3), atlas_roi.B));

annot.atlas_roi = annot.colortable.table(temp1,5);

roiVert = allvertices;
roiVert(~ismember(annot.label, annot.atlas_roi),:) = nan;

annot.vertices_roi = annot.vertices(ismember(annot.label, annot.atlas_roi));

temp2 = and(and((ismember(allfaces(:,1),annot.vertices_roi)),...
    (ismember(allfaces(:,2),annot.vertices_roi))),...
    sum(ismember(allfaces(:,3),annot.vertices_roi)));

roiFace = allfaces;
roiFace(~temp2,:) = nan;



end