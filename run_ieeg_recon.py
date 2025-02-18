#!/usr/bin/env python3

import argparse
from pathlib import Path
from src.ieeg_recon import run_pipeline

def main():
    parser = argparse.ArgumentParser(
        description='Run IEEG Reconstruction Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('--t1', type=str, required=True,
                       help='Path to pre-implant T1 MRI')
    parser.add_argument('--ct', type=str, required=True,
                       help='Path to post-implant CT')
    parser.add_argument('--elec', type=str, required=True,
                       help='Path to electrode coordinates file')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory path')

    # Optional arguments
    parser.add_argument('--freesurfer-dir', type=str,
                       help='Path to FreeSurfer directory (optional)')
    parser.add_argument('--env-path', type=str, default='.env',
                       help='Path to environment file')
    parser.add_argument('--modules', type=str, default='1,2,3',
                       help='Modules to run (comma-separated, e.g., "1,2,3")')
    parser.add_argument('--skip-existing', action='store_true',
                       help='Skip processing if output files exist')
    parser.add_argument('--reg-type', type=str, default='gc_noCTthereshold',
                       choices=['gc', 'g', 'gc_noCTthereshold'],
                       help='Registration type')
    parser.add_argument('--qa-viewer', type=str, default='niplot',
                       choices=['freeview', 'freeview_snapshot', 'niplot', 'itksnap', 'none'],
                       help='Quality assurance viewer type')

    args = parser.parse_args()

    # Convert paths to Path objects
    project_path = Path(__file__).parent
    pre_implant_mri = Path(args.t1)
    post_implant_ct = Path(args.ct)
    ct_electrodes = Path(args.elec)
    output_dir = Path(args.output_dir)
    env_path = project_path / args.env_path
    freesurfer_dir = Path(args.freesurfer_dir) if args.freesurfer_dir else None

    # Validate input files exist
    for filepath in [pre_implant_mri, post_implant_ct, ct_electrodes]:
        if not filepath.exists():
            print(f"Error: File not found: {filepath}")
            return

    # Convert modules string to list
    modules = args.modules.split(',')

    print("Processing with following inputs:")
    print(f"T1 MRI: {pre_implant_mri}")
    print(f"CT: {post_implant_ct}")
    print(f"Electrodes: {ct_electrodes}")
    print(f"Output directory: {output_dir}")
    print(f"FreeSurfer directory: {freesurfer_dir}")
    print(f"Using modules: {modules}")

    try:
        # Run pipeline with specified parameters
        run_pipeline(
            pre_implant_mri=pre_implant_mri,
            post_implant_ct=post_implant_ct,
            ct_electrodes=ct_electrodes,
            output_dir=output_dir,
            env_path=env_path,
            freesurfer_dir=freesurfer_dir,
            modules=modules,
            skip_existing=args.skip_existing,
            reg_type=args.reg_type,
            qa_viewer=args.qa_viewer
        )
        print("Processing complete!")
    except Exception as e:
        print(f"Error during pipeline execution: {str(e)}")

if __name__ == "__main__":
    main()