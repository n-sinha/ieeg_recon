#!/usr/bin/env python3
"""
ITK-SNAP Workspace Creator for iEEG Reconstruction

This script creates an ITK-SNAP workspace file (.itksnap) for visualizing
electrode spheres in the iEEG reconstruction pipeline.

Usage:
    python create_itksnap_workspace.py <output_dir> <pre_implant_mri> <ct_to_mri> <electrodes_inMRI> <electrode_names>
"""

import sys
import numpy as np
from pathlib import Path


def create_itksnap_workspace(output_dir, pre_implant_mri, ct_to_mri, electrodes_inMRI, electrode_names_file):
    """
    Create an ITK-SNAP workspace file for visualizing electrode spheres
    
    This function creates an ITK-SNAP workspace file (.itksnap) that includes:
    - Layer 0: Pre-implant MRI as the main anatomical image
    - Layer 1: Registered CT as an overlay (50% transparency)
    - Layer 2: Electrode spheres as segmentation labels
    
    Each electrode is assigned a unique color and label in the workspace.
    The workspace file can be opened directly in ITK-SNAP for interactive
    visualization of the electrode placement.
    
    Args:
        output_dir (str/Path): Output directory for module 2
        pre_implant_mri (str/Path): Path to pre-implant MRI
        ct_to_mri (str/Path): Path to registered CT image
        electrodes_inMRI (str/Path): Path to electrode spheres image
        electrode_names_file (str/Path): Path to electrode names file
    """
    # Convert to Path objects
    output_dir = Path(output_dir)
    pre_implant_mri = Path(pre_implant_mri)
    ct_to_mri = Path(ct_to_mri)
    electrodes_inMRI = Path(electrodes_inMRI)
    electrode_names_file = Path(electrode_names_file)
    
    # Load electrode names
    electrode_names = np.loadtxt(electrode_names_file, dtype=str)
    
    # Define file paths with quotes (as required by ITK-SNAP)
    mri_path = f'"{pre_implant_mri}"'
    ct_path = f'"{ct_to_mri}"'
    spheres_path = f'"{electrodes_inMRI}"'
    module2_quotes = f'"{output_dir}"'
    
    # Create the ITK-SNAP workspace XML content
    html_top = f'''<?xml version="1.0" encoding="UTF-8" ?>
<!--ITK-SNAP (itksnap.org) Project File

This file can be moved/copied along with the images that it references
as long as the relative location of the images to the project file is 
the same. Do not modify the SaveLocation entry, or this will not work.
-->
<!DOCTYPE registry [
<!ELEMENT registry (entry*,folder*)>
<!ELEMENT folder (entry*,folder*)>
<!ELEMENT entry EMPTY>
<!ATTLIST folder key CDATA #REQUIRED>
<!ATTLIST entry key CDATA #REQUIRED>
<!ATTLIST entry value CDATA #REQUIRED>
]>
<registry>
  <entry key="SaveLocation" value={module2_quotes} />
  <entry key="Version" value="20190612" />
  <folder key="Annotations" >
    <entry key="Format" value="ITK-SNAP Annotation File" />
    <entry key="FormatDate" value="20150624" />
  </folder>
  <folder key="Layers" >
    <folder key="Layer[000]" >
      <entry key="AbsolutePath" value={mri_path} />
      <entry key="Role" value="MainRole" />
      <entry key="Tags" value="" />
      <folder key="IOHints" >
      </folder>
      <folder key="LayerMetaData" >
        <entry key="Alpha" value="255" />
        <entry key="CustomNickName" value="" />
        <entry key="Sticky" value="0" />
        <entry key="Tags" value="" />
        <folder key="DisplayMapping" >
          <folder key="ColorMap" >
            <entry key="Preset" value="Grayscale" />
          </folder>
          <folder key="Curve" >
            <entry key="NumberOfControlPoints" value="3" />
            <folder key="ControlPoint[0]" >
              <entry key="tValue" value="0" />
              <entry key="xValue" value="0" />
            </folder>
            <folder key="ControlPoint[1]" >
              <entry key="tValue" value="0.5" />
              <entry key="xValue" value="0.5" />
            </folder>
            <folder key="ControlPoint[2]" >
              <entry key="tValue" value="1" />
              <entry key="xValue" value="1" />
            </folder>
          </folder>
        </folder>
      </folder>
      <folder key="ProjectMetaData" >
        <entry key="GaussianBlurScale" value="1" />
        <entry key="RemappingExponent" value="3" />
        <entry key="RemappingSteepness" value="0.04" />
        <folder key="Files" >
          <folder key="Grey" >
            <entry key="Dimensions" value="192 256 160" />
            <entry key="Orientation" value="LPI" />
          </folder>
        </folder>
        <folder key="IOHistory" >
          <folder key="AnatomicImage" >
            <entry key="ArraySize" value="1" />
            <entry key="Element[0]" value={ct_path} />
          </folder>
          <folder key="LabelImage" >
            <entry key="ArraySize" value="1" />
            <entry key="Element[0]" value={spheres_path} />
          </folder>
        </folder>
        <folder key="IRIS" >
          <entry key="SliceViewLayerLayout" value="Stacked" />
          <folder key="BoundingBox" >
            <entry key="InterpolationMethod" value="Nearest" />
            <entry key="ResampleDimensions" value="192 256 160" />
            <entry key="SeedWithCurrentSegmentation" value="0" />
            <folder key="ROIBox[0]" >
              <entry key="Index" value="0" />
              <entry key="Size" value="192" />
            </folder>
            <folder key="ROIBox[1]" >
              <entry key="Index" value="0" />
              <entry key="Size" value="256" />
            </folder>
            <folder key="ROIBox[2]" >
              <entry key="Index" value="0" />
              <entry key="Size" value="160" />
            </folder>
          </folder>
          <folder key="DisplayMapping" >
            <folder key="ColorMap" >
              <entry key="Preset" value="Grayscale" />
            </folder>
            <folder key="Curve" >
              <entry key="NumberOfControlPoints" value="3" />
              <folder key="ControlPoint[0]" >
                <entry key="tValue" value="0" />
                <entry key="xValue" value="0" />
              </folder>
              <folder key="ControlPoint[1]" >
                <entry key="tValue" value="0.5" />
                <entry key="xValue" value="0.5" />
              </folder>
              <folder key="ControlPoint[2]" >
                <entry key="tValue" value="1" />
                <entry key="xValue" value="1" />
              </folder>
            </folder>
          </folder>
          <folder key="LabelState" >
            <entry key="CoverageMode" value="OverAll" />
            <entry key="DrawingLabel" value="1" />
            <entry key="OverwriteLabel" value="0" />
            <entry key="PolygonInvert" value="0" />
            <entry key="SegmentationAlpha" value="0.5" />
          </folder>'''

    html_middle = f'''<folder key="LabelTable" >
            <entry key="NumberOfElements" value="{len(electrode_names)}" />
            '''

    for i, label_name in enumerate(electrode_names):
        html_middle_template = f'''<folder key="Element[{i}]" >
                <entry key="Alpha" value="255" />
                <entry key="Color" value="0 255 0" />
                <entry key="Flags" value="1 1" />
                <entry key="Index" value="{i+1}" />
                <entry key="Label" value="{label_name}" />
              </folder>'''
        html_middle = html_middle + html_middle_template

    html_bottom = f'''
</folder>
          <folder key="MeshOptions" >
            <entry key="DecimateFeatureAngle" value="45" />
            <entry key="DecimateMaximumError" value="0.002" />
            <entry key="DecimatePreserveTopology" value="1" />
            <entry key="DecimateTargetReduction" value="0.95" />
            <entry key="GaussianError" value="0.03" />
            <entry key="GaussianStandardDeviation" value="0.8" />
            <entry key="MeshSmoothingBoundarySmoothing" value="0" />
            <entry key="MeshSmoothingConvergence" value="0" />
            <entry key="MeshSmoothingFeatureAngle" value="45" />
            <entry key="MeshSmoothingFeatureEdgeSmoothing" value="0" />
            <entry key="MeshSmoothingIterations" value="20" />
            <entry key="MeshSmoothingRelaxationFactor" value="0.01" />
            <entry key="UseDecimation" value="0" />
            <entry key="UseGaussianSmoothing" value="1" />
            <entry key="UseMeshSmoothing" value="0" />
          </folder>
        </folder>
        <folder key="SNAP" >
          <folder key="SnakeParameters" >
            <entry key="AdvectionSpeedExponent" value="0" />
            <entry key="AdvectionWeight" value="0" />
            <entry key="AutomaticTimeStep" value="1" />
            <entry key="Clamp" value="1" />
            <entry key="CurvatureSpeedExponent" value="-1" />
            <entry key="CurvatureWeight" value="0.2" />
            <entry key="Ground" value="5" />
            <entry key="LaplacianSpeedExponent" value="0" />
            <entry key="LaplacianWeight" value="0" />
            <entry key="PropagationSpeedExponent" value="1" />
            <entry key="PropagationWeight" value="1" />
            <entry key="SnakeType" value="RegionCompetition" />
            <entry key="SolverAlgorithm" value="ParallelSparseField" />
            <entry key="TimeStepFactor" value="1" />
          </folder>
        </folder>
      </folder>
    </folder>
    <folder key="Layer[001]" >
      <entry key="AbsolutePath" value={ct_path} />
      <entry key="Role" value="OverlayRole" />
      <entry key="Tags" value="" />
      <folder key="IOHints" >
      </folder>
      <folder key="ImageTransform" >
        <entry key="IsIdentity" value="1" />
      </folder>
      <folder key="LayerMetaData" >
        <entry key="Alpha" value="0.5" />
        <entry key="CustomNickName" value="" />
        <entry key="Sticky" value="1" />
        <entry key="Tags" value="" />
        <folder key="DisplayMapping" >
          <folder key="ColorMap" >
            <entry key="Preset" value="Grayscale" />
          </folder>
          <folder key="Curve" >
            <entry key="NumberOfControlPoints" value="3" />
            <folder key="ControlPoint[0]" >
              <entry key="tValue" value="0" />
              <entry key="xValue" value="0" />
            </folder>
            <folder key="ControlPoint[1]" >
              <entry key="tValue" value="0.3125" />
              <entry key="xValue" value="0.5" />
            </folder>
            <folder key="ControlPoint[2]" >
              <entry key="tValue" value="0.625" />
              <entry key="xValue" value="1" />
            </folder>
          </folder>
        </folder>
      </folder>
    </folder>
    <folder key="Layer[002]" >
      <entry key="AbsolutePath" value={spheres_path} />
      <entry key="Role" value="SegmentationRole" />
      <entry key="Tags" value="" />
      <folder key="IOHints" >
      </folder>
      <folder key="LayerMetaData" >
        <entry key="Alpha" value="0" />
        <entry key="CustomNickName" value="" />
        <entry key="Sticky" value="1" />
        <entry key="Tags" value="" />
      </folder>
    </folder>
  </folder>
</registry>
'''

    html_all = html_top + html_middle + html_bottom

    # Save the workspace file
    workspace_file = output_dir / 'electrode_workspace.itksnap'
    with open(workspace_file, 'w') as f:
        f.write(html_all)
    
    print(f"ITK-SNAP workspace created: {workspace_file}")
    return str(workspace_file)


def main():
    """Main function for command-line usage"""
    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(1)
    
    output_dir = sys.argv[1]
    pre_implant_mri = sys.argv[2]
    ct_to_mri = sys.argv[3]
    electrodes_inMRI = sys.argv[4]
    electrode_names_file = sys.argv[5]
    
    try:
        workspace_file = create_itksnap_workspace(
            output_dir, pre_implant_mri, ct_to_mri, electrodes_inMRI, electrode_names_file
        )
        print(f"Successfully created workspace: {workspace_file}")
    except Exception as e:
        print(f"Error creating workspace: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 