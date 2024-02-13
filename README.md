# iEEG_implant_reconstruction

## Startup 
* Make sure that FSL is setup properly and added to the shell profile.

## Download example data

* To test the codes, download example data and put it in this directory: [Dropbox](https://www.dropbox.com/sh/ylxc586grm0p7au/AAAs8QQwUo0VQOSweDyj1v_ta?dl=0)


## Module 1

"module1" exports the electrode coordinates of a post-implant CT scan in both voxel and native spaces. The function uses the information from a config file "config_iEEGrecon" and reads in the electrode coordinates from a file specified by "obj.postImplantCT_electrodes". The function then creates a directory "ieeg_recon/module1" under the "obj.output" directory. The electrode names and their corresponding coordinates in CT voxel space are saved in two separate text files, "electrode_names.txt" and "electrodes_inCTvox.txt". The CT image data and its header information are read using the "niftiinfo" and "niftiread" functions and stored in the "CT" structure. The function removes the skull portion of the CT image by finding voxels with intensity values higher than a threshold (99.95 percentile). The electrode coordinates in CT voxel space and CT native space (in mm) are calculated and saved in two separate text files, "electrodes_inCTmm.txt".

## Module 2
This is a function in MATLAB (or similar) programming language. The purpose of this function is to perform image registration between two images: the pre-implant MRI and post-implant CT scans. The function is part of a larger program called iEEGrecon. The function makes use of external tools such as FSL and ITK-SNAP to perform the image registration. The function first checks if the required output files are present, and if not, performs the image registration. The output files include the registered CT image, transformation matrix, and electrode locations in both MRI and CT images. The function saves the output files in a folder named ieeg_recon/module2 within the output folder. The function also has a try-catch block to handle potential errors that might occur during the image registration.

This code segment is performing image registration between two medical images: a post-implant CT scan and a pre-implant MRI scan. The code uses two registration methods: "greedy with image centering" and "FLIRT registration fine-tuned with greedy".

In the first method, "greedy with image centering", the CT scan is thresholded with a value of 0, then registered to the MRI using the 'greedy' program with "image centering" option enabled. The output of this registration is a transform matrix that is converted to a text file in FSL format. This transform is then applied to the thresholded CT image to align it with the MRI.

In the second method, "FLIRT registration fine-tuned with greedy", the thresholded CT is first registered to the MRI using the FSL's 'flirt' program with mutual information cost function and 6 degrees of freedom. The output of this registration is then fine-tuned using the 'greedy' program with "identity" option enabled.

Both methods produce an aligned CT image in the end, which can be used for further processing or analysis.

## Module 2 Quality Assurance
The code is a MATLAB function for quality assurance (QA) in module 2 of a iEEG reconstruction process. It creates a folder "ieeg_recon/module2" in the specified "obj.output" directory and generates screenshots of the registration of the CT and MRI images using the specified image viewer. The screenshots include views of the registration in sagittal, coronal, axial, and 3D views. The image viewer can be specified as "freeview", "freeview_snapshot", or "itksnap". If the image viewer is not specified, the code throws an error. The screenshots are saved in the "ieeg_recon/module2" folder.

## Module 3
This code defines a function module3 that maps electrodes in a patient's brain to regions of interest (ROIs) defined by an atlas and a lookup table.

The function takes three inputs: obj, atlas, and lookupTable. obj is an object with a property output that specifies the output directory. atlas is the file path to an atlas image in NIfTI format. lookupTable is a table that maps values in the atlas image to ROI labels.

The code first creates the output directory if it does not already exist. Then, it loads electrode coordinates in both millimeters and voxels from text files in the obj.output directory. The atlas and lookup table are also loaded into the workspace.

The code converts the voxel coordinates of the atlas to millimeter coordinates using the transformation matrix in the NIfTI header. Then, it maps each electrode to the closest ROI using a k-nearest neighbors search with a distance threshold of 2.6 millimeters. Electrodes that are farther than 2.6 millimeters from any ROI are considered to be in the white matter or outside the brain.

The final result is a table electrodes2ROI that contains the label, coordinates, and ROI information for each electrode. The table is saved to a .csv file in the output directory.

## Outline
Here is an outline for the README documentation of your code:

### Introduction:

Brief overview of the purpose and functionality of the code.
Description of what the code does and its intended use.
Requirements:

List of required software and libraries (e.g. MATLAB, Python libraries).
Description of system requirements (e.g. operating system, hardware specifications).
Installation:

Step-by-step instructions for installing and setting up the code.
Explanation of how to set up dependencies and required libraries.
Usage:

Description of how to run the code and its input parameters.
Explanation of how to use the code and how to interpret its outputs.
Code Structure:

Explanation of the code's structure and organization.
List of files included in the code and their functions.
Explanation of the main functions and subroutines of the code.
Data Inputs:

Description of the input data format and requirements.
Explanation of how to prepare input data and where to obtain it.
Outputs:

Description of the code's outputs and the format of the output files.
Explanation of how to interpret the outputs and what they represent.
Limitations:

Explanation of any known limitations or problems with the code.
Discussion of any potential issues or limitations that users should be aware of.
References:

List of any relevant references or resources used in creating the code.
Explanation of how the code builds upon or relates to existing research or projects.
Contact:

Information on how to contact the authors or developers of the code.
Description of how to report bugs, ask questions, or provide feedback.

License:
Information on the license under which the code is released.
Explanation of the terms and conditions of using the code.

Contributing:
Information on how to contribute to the development of the code.
Explanation of the guidelines for contributing and how to get involved.