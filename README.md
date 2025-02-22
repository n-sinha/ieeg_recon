# iEEG Implant Reconstruction Pipeline

A Python-based pipeline for reconstructing intracranial EEG (iEEG) electrode locations from post-implant CT and pre-implant MRI scans.

## Overview

This pipeline provides tools for:
- Extracting electrode coordinates from post-implant CT scans
- Co-registering CT and MRI images
- Mapping electrodes to anatomical regions
- Quality assurance visualization

## Prerequisites

### Required Software
- FSL (FMRIB Software Library)
- ITK-SNAP
- FreeSurfer
- Python 3.x with required packages (see requirements.txt)

### System Requirements
- Operating System: Linux/Unix (recommended), macOS, or Windows with Unix subsystem
- RAM: 8GB minimum, 16GB recommended
- Storage: At least 10GB free space for processing

## Installation

There are two ways to install this package:

1. **Direct installation from GitHub**:
   ```bash
   pip install git+https://github.com/n-sinha/ieeg_recon.git
   ```
   
   For development installation (includes testing tools):
   ```bash
   pip install "git+https://github.com/n-sinha/ieeg_recon.git#egg=ieeg_recon[dev]"
   ```

2. **Local development installation**:
   ```bash
   # Clone this repository
   git clone https://github.com/n-sinha/ieeg_recon.git
   cd ieeg_recon
   
   # Install in editable mode with development dependencies
   pip install -e ".[dev]"
   ```

3. Create and configure your `.env` file:
   ```bash
   # FSL Configuration
   FSL_DIR=/path/to/fsl
   FSL_OUTPUT_TYPE=NIFTI_GZ

   # ITK-SNAP Configuration
   ITKSNAP_DIR=/path/to/ITK-SNAP

   # FreeSurfer Configuration
   FREESURFER_HOME=/path/to/freesurfer
   SUBJECTS_DIR=/path/to/subjects_dir
   ```

   Example of a configured .env file:
   ```bash
   # FSL Configuration
   FSL_DIR=/usr/local/fsl
   FSL_OUTPUT_TYPE=NIFTI_GZ

   # ITK-SNAP Configuration
   ITKSNAP_DIR=/Applications/ITK-SNAP.app/Contents/bin

   # FreeSurfer Configuration
   FREESURFER_HOME=/Applications/freesurfer/7.3.2
   SUBJECTS_DIR=/Applications/freesurfer/7.3.2/subjects
   ```

After installation, you can import the package in Python:
```python
from ieeg_recon import IEEGRecon, run_pipeline
```
## Usage

The pipeline can be run using the command-line interface:

```bash
python run_ieeg_recon.py --t1 /path/to/t1.nii.gz \
                        --ct /path/to/ct.nii.gz \
                        --elec /path/to/electrodes.txt \
                        --output-dir /path/to/output \
                        --freesurfer-dir /path/to/freesurfer_dir
```

### Command Line Arguments

Required:
- `--t1`: Path to pre-implant T1 MRI
- `--ct`: Path to post-implant CT
- `--elec`: Path to electrode coordinates file
- `--output-dir`: Output directory path

Optional:
- `--freesurfer-dir`: Path to FreeSurfer directory
- `--env-path`: Path to environment file (default: '.env')
- `--modules`: Modules to run (comma-separated, e.g., "1,2,3,4", default: "1,2,3,4")
- `--skip-existing`: Skip processing if output files exist
- `--reg-type`: Registration type ('gc', 'g', 'gc_noCTthereshold', default: 'gc_noCTthereshold')
- `--qa-viewer`: Quality assurance viewer type ('freeview', 'freeview_snapshot', 'niplot', 'itksnap', 'none', default: 'niplot')

### Pipeline Modules

#### Module 1: Electrode Coordinate Export
Exports electrode coordinates from post-implant CT in voxel and native space.

#### Module 2: Image Registration
Performs CT-MRI co-registration with options for:
- Greedy registration with image centering ('gc')
- FLIRT registration fine-tuned with greedy ('g')
- Greedy registration without CT thresholding ('gc_noCTthereshold')

#### Module 3: ROI Mapping
Maps electrodes to anatomical regions using an atlas and provides quality assurance visualization.

#### Module 4: MNI Space Transformation
Transforms electrode coordinates to MNI305 and MNI152 standard spaces for cross-subject comparisons and group analyses.

## Output Structure

```
output_dir/
├── module1/
│   ├── electrode_names.txt       # List of electrode names
│   ├── electrodes_inCTvox.txt   # Electrode coordinates in CT voxel space
│   └── electrodes_inCTmm.txt    # Electrode coordinates in CT millimeter space
├── module2/
│   ├── ct_to_mri.nii.gz         # Registered CT to MRI image
│   ├── electrodes_inMRI.nii.gz  # Electrode positions in MRI space
│   └── QA_registration.png      # Quality assurance visualization
└── module3/
    ├── electrodes2ROI.csv       # Mapping of electrodes to regions of interest
    └── electrode_visualization.html  # Interactive 3D visualization
```
## Quality Assurance

The pipeline includes comprehensive quality assurance tools:
- Interactive 3D visualization of electrode positions
- Multiple viewing options (hemisphere visibility, electrode size, label toggle)
- Registration quality checks
- Electrode placement verification

## Citation

If you use this software in your research, please cite:

> **Lucas A, Scheid BH, Pattnaik AR, Gallagher R, Mojena M, Tranquille A, Prager B, Gleichgerrcht E, Gong R, Litt B, Davis KA, Das S, Stein JM, Sinha N. iEEG-recon: A fast and scalable pipeline for accurate reconstruction of intracranial electrodes and implantable devices. Epilepsia. 2024 Mar;65(3):817-829. doi: 10.1111/epi.17863.**

## Support

For bug reports and feature requests, please use the [GitHub Issue Tracker](https://github.com/n-sinha/ieeg_recon/issues).


