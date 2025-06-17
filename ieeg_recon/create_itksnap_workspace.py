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


def create_itksnap_label_file(electrode_names, output_file):
    """
    Create an ITK-SNAP label description file with individual colors for each electrode
    
    Args:
        electrode_names (list): List of electrode names
        output_file (Path): Output file path for the label file
    """
    # Define a color palette for electrodes (similar to the reference file)
    colors = [
        (30, 104, 225),   # Blue
        (218, 37, 211),   # Magenta
        (171, 94, 31),    # Brown
        (156, 122, 3),    # Olive
        (133, 3, 78),     # Purple
        (65, 186, 154),   # Teal
        (207, 245, 181),  # Light green
        (189, 42, 141),   # Pink
        (201, 188, 30),   # Yellow
        (136, 224, 245),  # Light blue
        (103, 182, 150),  # Mint
        (17, 210, 152),   # Green
        (57, 126, 27),    # Dark green
        (58, 5, 32),      # Dark red
        (135, 29, 9),     # Red
        (189, 247, 183),  # Light mint
        (72, 179, 247),   # Sky blue
        (165, 14, 226),   # Violet
        (229, 194, 54),   # Gold
        (153, 132, 25),   # Dark yellow
        (40, 165, 20),    # Forest green
        (145, 226, 11),   # Lime
        (153, 174, 77),   # Olive green
        (116, 67, 94),    # Plum
        (239, 151, 92),   # Orange
        (31, 233, 24),    # Bright green
        (219, 114, 246),  # Light purple
        (232, 8, 136),    # Hot pink
        (182, 47, 116),   # Rose
        (155, 153, 215),  # Lavender
        (72, 102, 51),    # Dark olive
        (1, 127, 254),    # Blue
        (122, 167, 52),   # Green
        (250, 220, 66),   # Yellow
        (4, 1, 128),      # Navy
        (227, 146, 169),  # Pink
        (176, 254, 167),  # Light green
        (32, 64, 86),     # Dark blue
        (94, 67, 15),     # Brown
        (247, 97, 65),    # Coral
        (245, 120, 129),  # Salmon
        (110, 177, 232),  # Light blue
        (10, 250, 73),    # Bright green
        (116, 28, 169),   # Purple
        (246, 29, 225),   # Magenta
        (117, 185, 53),   # Green
        (169, 130, 168),  # Lavender
        (180, 165, 174),  # Gray
        (79, 135, 212),   # Blue
        (99, 105, 45),    # Olive
        (254, 254, 205),  # Cream
        (99, 91, 205),    # Periwinkle
        (179, 75, 41),    # Brown
        (16, 235, 109),   # Green
        (70, 98, 200),    # Blue
        (134, 71, 142),   # Purple
        (209, 195, 186),  # Beige
        (243, 156, 4),    # Orange
        (67, 72, 17),     # Dark green
        (150, 194, 30),   # Green
        (98, 186, 167),   # Teal
        (94, 59, 15),     # Brown
        (36, 118, 159),   # Blue
        (220, 153, 86),   # Tan
        (15, 202, 222),   # Cyan
        (197, 113, 16),   # Brown
        (23, 37, 30),     # Dark green
        (68, 115, 100),   # Green
        (84, 144, 222),   # Blue
        (106, 71, 235),   # Purple
        (59, 64, 185),    # Blue
        (186, 113, 47),   # Brown
        (245, 164, 162),  # Pink
        (3, 147, 239),    # Blue
        (101, 177, 23),   # Green
        (168, 55, 167),   # Purple
        (176, 235, 238),  # Light blue
        (147, 18, 49),    # Red
        (50, 57, 179),    # Blue
        (176, 238, 207),  # Mint
        (131, 114, 239),  # Purple
        (64, 5, 145),     # Purple
        (186, 44, 228),   # Magenta
        (155, 161, 85),   # Olive
        (81, 169, 20),    # Green
        (225, 211, 35),   # Yellow
        (163, 51, 57),    # Red
        (97, 215, 1),     # Green
        (156, 240, 152),  # Light green
        (185, 224, 108),  # Light green
        (218, 216, 174),  # Beige
        (129, 166, 53),   # Green
        (100, 167, 106),  # Green
        (165, 177, 103),  # Green
        (58, 160, 25),    # Green
        (15, 3, 111),     # Blue
        (19, 237, 198),   # Cyan
        (130, 13, 131),   # Purple
        (139, 250, 231),  # Light cyan
        (15, 182, 230),   # Blue
        (202, 4, 109),    # Purple
        (136, 63, 198),   # Purple
        (190, 44, 194),   # Pink
        (40, 112, 79),    # Green
        (19, 125, 102),   # Green
        (137, 162, 84),   # Olive
        (78, 207, 52),    # Green
        (95, 33, 218),    # Purple
        (250, 233, 215),  # Cream
        (225, 158, 84),   # Tan
        (189, 21, 45),    # Red
        (40, 250, 144),   # Green
        (219, 249, 143),  # Light green
        (4, 45, 228),     # Blue
        (96, 220, 57),    # Green
        (130, 34, 192),   # Purple
        (35, 156, 224),   # Blue
        (154, 66, 139),   # Purple
        (4, 56, 161),     # Blue
        (73, 103, 226),   # Blue
        (86, 137, 158),   # Teal
        (32, 235, 35),    # Green
        (101, 127, 67),   # Olive
        (128, 22, 216),   # Purple
        (229, 84, 28),    # Orange
        (61, 244, 244),   # Cyan
        (71, 235, 113),   # Green
        (174, 209, 36),   # Green
        (37, 111, 30),    # Green
        (145, 223, 57),   # Green
        (42, 33, 194),    # Blue
        (74, 129, 28),    # Green
        (241, 254, 183),  # Light green
        (52, 129, 124),   # Teal
        (2, 227, 108),    # Green
        (248, 35, 36),    # Red
        (156, 104, 238),  # Purple
        (189, 5, 245),    # Magenta
        (191, 10, 228),   # Purple
        (55, 8, 187),     # Blue
        (94, 188, 81),    # Green
        (150, 234, 219),  # Light cyan
        (132, 33, 92),    # Purple
        (34, 79, 195),    # Blue
        (241, 66, 54),    # Red
        (72, 136, 181),   # Blue
        (52, 105, 205),   # Blue
        (10, 126, 220),   # Blue
        (120, 14, 223),   # Purple
        (168, 241, 13),   # Green
        (127, 229, 49),   # Green
        (70, 144, 179),   # Blue
        (92, 171, 190),   # Blue
        (85, 50, 143),    # Purple
        (79, 102, 90),    # Gray
        (177, 140, 193),  # Purple
        (217, 174, 121),  # Tan
        (65, 216, 157),   # Green
        (104, 76, 65),    # Brown
        (191, 1, 243),    # Magenta
        (224, 103, 163),  # Pink
        (102, 119, 34),   # Green
        (180, 54, 67),    # Red
        (36, 109, 168),   # Blue
        (235, 137, 216),  # Pink
        (46, 252, 181),   # Light green
        (40, 157, 84),    # Green
        (10, 219, 69),    # Green
        (120, 246, 53),   # Green
        (254, 118, 49),   # Orange
        (151, 124, 56),   # Brown
        (168, 124, 155),  # Pink
        (95, 40, 138),    # Purple
        (122, 254, 165),  # Light green
        (45, 120, 7),     # Green
        (132, 104, 37),   # Brown
        (236, 74, 159),   # Pink
        (118, 8, 87),     # Purple
        (63, 23, 82),     # Purple
        (145, 225, 184),  # Light green
        (204, 221, 187),  # Light gray
        (232, 8, 60),     # Red
        (59, 64, 201),    # Blue
        (145, 77, 2),     # Brown
        (83, 14, 242),    # Blue
        (237, 3, 203),    # Magenta
        (60, 152, 27),    # Green
        (149, 43, 219),   # Purple
        (239, 7, 35),     # Red
        (88, 131, 115),   # Gray
        (130, 55, 0),     # Brown
        (254, 139, 106),  # Orange
        (206, 144, 186),  # Pink
        (220, 5, 178),    # Magenta
        (201, 228, 122),  # Light green
        (118, 143, 134),  # Gray
        (59, 70, 16),     # Dark green
        (135, 182, 91),   # Green
        (9, 36, 187),     # Blue
        (92, 82, 218),    # Purple
        (44, 253, 122),   # Green
        (178, 101, 115),  # Pink
        (76, 65, 119),    # Purple
        (217, 61, 78),    # Red
        (124, 118, 56),   # Brown
        (41, 92, 13),     # Green
        (34, 143, 154),   # Teal
        (225, 75, 9),     # Orange
        (182, 156, 179),  # Gray
        (134, 101, 179),  # Purple
        (43, 189, 134),   # Green
        (27, 155, 226),   # Blue
        (5, 122, 123),    # Teal
        (61, 58, 101),    # Blue
        (158, 172, 141),  # Gray
        (24, 32, 196),    # Blue
        (119, 18, 239),   # Purple
        (79, 253, 23),    # Green
        (123, 179, 17),   # Green
        (142, 213, 144),  # Light green
        (199, 187, 136),  # Beige
        (73, 80, 169),    # Blue
        (20, 254, 208),   # Cyan
        (69, 230, 220),   # Cyan
    ]
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("################################################\n")
        f.write("# ITK-SnAP Label Description File\n")
        f.write("# File format: \n")
        f.write("# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL\n")
        f.write("# Fields: \n")
        f.write("#    IDX:   Zero-based index \n")
        f.write("#    -R-:   Red color component (0..255)\n")
        f.write("#    -G-:   Green color component (0..255)\n")
        f.write("#    -B-:   Blue color component (0..255)\n")
        f.write("#    -A-:   Label transparency (0.00 .. 1.00)\n")
        f.write("#    VIS:   Label visibility (0 or 1)\n")
        f.write("#    IDX:   Label mesh visibility (0 or 1)\n")
        f.write("#  LABEL:   Label description \n")
        f.write("################################################\n")
        
        # Write clear label (index 0)
        f.write('0\t0\t0\t0\t0\t0\t0\t"Clear Label"\n')
        
        # Write electrode labels
        for i, name in enumerate(electrode_names):
            if i < len(colors):
                r, g, b = colors[i]
            else:
                # If we run out of colors, cycle through them
                r, g, b = colors[i % len(colors)]
            
            f.write(f'{i+1}\t{r}\t{g}\t{b}\t1\t1\t1\t"{name}"\n')
    
    print(f"ITK-SNAP label file created: {output_file}")


def create_itksnap_workspace(output_dir, pre_implant_mri, ct_to_mri, electrodes_inMRI, electrode_names_file):
    """
    Create an ITK-SNAP workspace file for visualizing electrode spheres
    
    This function creates an ITK-SNAP workspace file (.itksnap) that includes:
    - Layer 0: Pre-implant MRI as the main anatomical image
    - Layer 1: Registered CT as an overlay (50% transparency)
    - Layer 2: Electrode spheres as segmentation labels with proper names
    
    Each electrode is assigned a unique color and proper label in the workspace.
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
    
    # Create ITK-SNAP label file with individual colors for each electrode
    label_file = output_dir / 'electrodes_itk_snap_labels.txt'
    create_itksnap_label_file(electrode_names, label_file)
    
    # Define file paths with quotes (as required by ITK-SNAP)
    mri_path = f'"{pre_implant_mri}"'
    ct_path = f'"{ct_to_mri}"'
    spheres_path = f'"{electrodes_inMRI}"'
    spheres_labels = f'"{label_file}"'
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
        <entry key="ShowLabels" value="1" />
        <entry key="LabelOpacity" value="1.0" />
        <entry key="LabelSize" value="12" />
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
          <folder key="LabelDescriptions" >
            <entry key="ArraySize" value="1" />
            <entry key="Element[0]" value={spheres_labels} />
          </folder>
          <folder key="LabelImage" >
            <entry key="ArraySize" value="1" />
            <entry key="Element[0]" value={spheres_path} />
          </folder>
        </folder>
        <folder key="IRIS" >
          <entry key="SliceViewLayerLayout" value="Stacked" />
          <entry key="LabelDescriptionsFile" value={spheres_labels} />
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
            <entry key="ShowLabels" value="1" />
            <entry key="LabelOpacity" value="1.0" />
            <entry key="LabelSize" value="12" />
          </folder>'''

    html_middle = f'''<folder key="LabelTable" >
            <entry key="NumberOfElements" value="{len(electrode_names)}" />
            '''

    # Define a color palette for electrodes (same as in create_itksnap_label_file)
    colors = [
        (30, 104, 225),   # Blue
        (218, 37, 211),   # Magenta
        (171, 94, 31),    # Brown
        (156, 122, 3),    # Olive
        (133, 3, 78),     # Purple
        (65, 186, 154),   # Teal
        (207, 245, 181),  # Light green
        (189, 42, 141),   # Pink
        (201, 188, 30),   # Yellow
        (136, 224, 245),  # Light blue
        (103, 182, 150),  # Mint
        (17, 210, 152),   # Green
        (57, 126, 27),    # Dark green
        (58, 5, 32),      # Dark red
        (135, 29, 9),     # Red
        (189, 247, 183),  # Light mint
        (72, 179, 247),   # Sky blue
        (165, 14, 226),   # Violet
        (229, 194, 54),   # Gold
        (153, 132, 25),   # Dark yellow
        (40, 165, 20),    # Forest green
        (145, 226, 11),   # Lime
        (153, 174, 77),   # Olive green
        (116, 67, 94),    # Plum
        (239, 151, 92),   # Orange
        (31, 233, 24),    # Bright green
        (219, 114, 246),  # Light purple
        (232, 8, 136),    # Hot pink
        (182, 47, 116),   # Rose
        (155, 153, 215),  # Lavender
        (72, 102, 51),    # Dark olive
        (1, 127, 254),    # Blue
        (122, 167, 52),   # Green
        (250, 220, 66),   # Yellow
        (4, 1, 128),      # Navy
        (227, 146, 169),  # Pink
        (176, 254, 167),  # Light green
        (32, 64, 86),     # Dark blue
        (94, 67, 15),     # Brown
        (247, 97, 65),    # Coral
        (245, 120, 129),  # Salmon
        (110, 177, 232),  # Light blue
        (10, 250, 73),    # Bright green
        (116, 28, 169),   # Purple
        (246, 29, 225),   # Magenta
        (117, 185, 53),   # Green
        (169, 130, 168),  # Lavender
        (180, 165, 174),  # Gray
        (79, 135, 212),   # Blue
        (99, 105, 45),    # Olive
        (254, 254, 205),  # Cream
        (99, 91, 205),    # Periwinkle
        (179, 75, 41),    # Brown
        (16, 235, 109),   # Green
        (70, 98, 200),    # Blue
        (134, 71, 142),   # Purple
        (209, 195, 186),  # Beige
        (243, 156, 4),    # Orange
        (67, 72, 17),     # Dark green
        (150, 194, 30),   # Green
        (98, 186, 167),   # Teal
        (94, 59, 15),     # Brown
        (36, 118, 159),   # Blue
        (220, 153, 86),   # Tan
        (15, 202, 222),   # Cyan
        (197, 113, 16),   # Brown
        (23, 37, 30),     # Dark green
        (68, 115, 100),   # Green
        (84, 144, 222),   # Blue
        (106, 71, 235),   # Purple
        (59, 64, 185),    # Blue
        (186, 113, 47),   # Brown
        (245, 164, 162),  # Pink
        (3, 147, 239),    # Blue
        (101, 177, 23),   # Green
        (168, 55, 167),   # Purple
        (176, 235, 238),  # Light blue
        (147, 18, 49),    # Red
        (50, 57, 179),    # Blue
        (176, 238, 207),  # Mint
        (131, 114, 239),  # Purple
        (64, 5, 145),     # Purple
        (186, 44, 228),   # Magenta
        (155, 161, 85),   # Olive
        (81, 169, 20),    # Green
        (225, 211, 35),   # Yellow
        (163, 51, 57),    # Red
        (97, 215, 1),     # Green
        (156, 240, 152),  # Light green
        (185, 224, 108),  # Light green
        (218, 216, 174),  # Beige
        (129, 166, 53),   # Green
        (100, 167, 106),  # Green
        (165, 177, 103),  # Green
        (58, 160, 25),    # Green
        (15, 3, 111),     # Blue
        (19, 237, 198),   # Cyan
        (130, 13, 131),   # Purple
        (139, 250, 231),  # Light cyan
        (15, 182, 230),   # Blue
        (202, 4, 109),    # Purple
        (136, 63, 198),   # Purple
        (190, 44, 194),   # Pink
        (40, 112, 79),    # Green
        (19, 125, 102),   # Green
        (137, 162, 84),   # Olive
        (78, 207, 52),    # Green
        (95, 33, 218),    # Purple
        (250, 233, 215),  # Cream
        (225, 158, 84),   # Tan
        (189, 21, 45),    # Red
        (40, 250, 144),   # Green
        (219, 249, 143),  # Light green
        (4, 45, 228),     # Blue
        (96, 220, 57),    # Green
        (130, 34, 192),   # Purple
        (35, 156, 224),   # Blue
        (154, 66, 139),   # Purple
        (4, 56, 161),     # Blue
        (73, 103, 226),   # Blue
        (86, 137, 158),   # Teal
        (32, 235, 35),    # Green
        (101, 127, 67),   # Olive
        (128, 22, 216),   # Purple
        (229, 84, 28),    # Orange
        (61, 244, 244),   # Cyan
        (71, 235, 113),   # Green
        (174, 209, 36),   # Green
        (37, 111, 30),    # Green
        (145, 223, 57),   # Green
        (42, 33, 194),    # Blue
        (74, 129, 28),    # Green
        (241, 254, 183),  # Light green
        (52, 129, 124),   # Teal
        (2, 227, 108),    # Green
        (248, 35, 36),    # Red
        (156, 104, 238),  # Purple
        (189, 5, 245),    # Magenta
        (191, 10, 228),   # Purple
        (55, 8, 187),     # Blue
        (94, 188, 81),    # Green
        (150, 234, 219),  # Light cyan
        (132, 33, 92),    # Purple
        (34, 79, 195),    # Blue
        (241, 66, 54),    # Red
        (72, 136, 181),   # Blue
        (52, 105, 205),   # Blue
        (10, 126, 220),   # Blue
        (120, 14, 223),   # Purple
        (168, 241, 13),   # Green
        (127, 229, 49),   # Green
        (70, 144, 179),   # Blue
        (92, 171, 190),   # Blue
        (85, 50, 143),    # Purple
        (79, 102, 90),    # Gray
        (177, 140, 193),  # Purple
        (217, 174, 121),  # Tan
        (65, 216, 157),   # Green
        (104, 76, 65),    # Brown
        (191, 1, 243),    # Magenta
        (224, 103, 163),  # Pink
        (102, 119, 34),   # Green
        (180, 54, 67),    # Red
        (36, 109, 168),   # Blue
        (235, 137, 216),  # Pink
        (46, 252, 181),   # Light green
        (40, 157, 84),    # Green
        (10, 219, 69),    # Green
        (120, 246, 53),   # Green
        (254, 118, 49),   # Orange
        (151, 124, 56),   # Brown
        (168, 124, 155),  # Pink
        (95, 40, 138),    # Purple
        (122, 254, 165),  # Light green
        (45, 120, 7),     # Green
        (132, 104, 37),   # Brown
        (236, 74, 159),   # Pink
        (118, 8, 87),     # Purple
        (63, 23, 82),     # Purple
        (145, 225, 184),  # Light green
        (204, 221, 187),  # Light gray
        (232, 8, 60),     # Red
        (59, 64, 201),    # Blue
        (145, 77, 2),     # Brown
        (83, 14, 242),    # Blue
        (237, 3, 203),    # Magenta
        (60, 152, 27),    # Green
        (149, 43, 219),   # Purple
        (239, 7, 35),     # Red
        (88, 131, 115),   # Gray
        (130, 55, 0),     # Brown
        (254, 139, 106),  # Orange
        (206, 144, 186),  # Pink
        (220, 5, 178),    # Magenta
        (201, 228, 122),  # Light green
        (118, 143, 134),  # Gray
        (59, 70, 16),     # Dark green
        (135, 182, 91),   # Green
        (9, 36, 187),     # Blue
        (92, 82, 218),    # Purple
        (44, 253, 122),   # Green
        (178, 101, 115),  # Pink
        (76, 65, 119),    # Purple
        (217, 61, 78),    # Red
        (124, 118, 56),   # Brown
        (41, 92, 13),     # Green
        (34, 143, 154),   # Teal
        (225, 75, 9),     # Orange
        (182, 156, 179),  # Gray
        (134, 101, 179),  # Purple
        (43, 189, 134),   # Green
        (27, 155, 226),   # Blue
        (5, 122, 123),    # Teal
        (61, 58, 101),    # Blue
        (158, 172, 141),  # Gray
        (24, 32, 196),    # Blue
        (119, 18, 239),   # Purple
        (79, 253, 23),    # Green
        (123, 179, 17),   # Green
        (142, 213, 144),  # Light green
        (199, 187, 136),  # Beige
        (73, 80, 169),    # Blue
        (20, 254, 208),   # Cyan
        (69, 230, 220),   # Cyan
    ]

    for i, label_name in enumerate(electrode_names):
        if i < len(colors):
            r, g, b = colors[i]
        else:
            # If we run out of colors, cycle through them
            r, g, b = colors[i % len(colors)]
        
        html_middle_template = f'''<folder key="Element[{i}]" >
                <entry key="Alpha" value="255" />
                <entry key="Color" value="{r} {g} {b}" />
                <entry key="Flags" value="1 1" />
                <entry key="Index" value="{i+1}" />
                <entry key="Label" value="{label_name}" />
                <entry key="Visible" value="1" />
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
        <entry key="ShowLabels" value="1" />
        <entry key="LabelOpacity" value="1.0" />
        <entry key="LabelSize" value="12" />
        <entry key="LabelDescriptionsFile" value={spheres_labels} />
        <entry key="AutoLoadLabels" value="1" />
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