#!/usr/bin/env python3

# Change log:
# v1.2.0. 2025-06-03. No changes, just tracking version with cellsnp_vcf_merger_submit.sh. 
# v1.1.0. 2025-06-03. No changes, just tracking version with cellsnp_vcf_merger_submit.sh. 
# v1.0.12. 2025-06-03. No changes, just tracking version with cellsnp_vcf_merger_submit.sh. 
# v1.0.11. 2025-06-03. No changes, just tracking version with cellsnp_vcf_merger_submit.sh. 
# v1.0.10. 2025-06-03. No changes, just tracking version with cellsnp_vcf_merger_submit.sh. 
# v1.0.9. 2025-06-03.
# v1.0.8. 2025-06-03.
# v1.0.7. 2025-06-02. No changes, just tracking version with cellsnp_vcf_merger_submit.sh.
# v1.0.6. 2025-06-02.
# v1.0.5. 2025-06-02.
# v1.0.4. 2025-06-02.
# v1.0.3. 2025-06-02. No changes, just tracking version with cellsnp_vcf_merger_submit.sh.
# v1.0.2. 2025-06-02. No changes, just tracking version with cellsnp_vcf_merger_submit.sh.
# v1.0.1. 2025-06-02.
# v1.0.0. 2025-04-01. Initial version.

# Version information
VERSION_NAME = 'CellSNP VCF Merger'
VERSION = '1.2.0'
VERSION_DATE = '2025-06-03'
COPYRIGHT = 'Copyright 1979-2025. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
COPYRIGHT_TEXT = '''
The MIT License (MIT).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
'''
DESCRIPTION = f'Merges all *_cellSNP_sorted.cells.vcf.gz files from a given study directory, then sorts and indexes the merged VCF using bcftools.\n\n{COPYRIGHT_TEXT}'

# --- Import necessary libraries ---
# General imports
import os
import sys
import subprocess
import logging
from pathlib import Path
import importlib

# Time and date
import time
from datetime import datetime, timedelta

# Argument parsing
import argparse

# --- Function to setup the logger ---
def setup_logger(logfilename, verbose=False, debug=False, output_dir=None):
    """
    Setup the logger to log to a file with the date and script name, and also log to the console.

    Arguments:
        logfilename (str): The base name for the log file.
        verbose (bool): If True, set logging level to DEBUG; otherwise, set to INFO.
        debug (bool): If True, set logging level to DEBUG; otherwise, set to INFO.
        output_dir (str): Directory where the log file will be saved. If None, uses current directory.
    Returns:
        logging.Logger: Configured logger instance.
    """

    # Get current date in the format YYYYMMDD
    date_str = datetime.now().strftime('%Y%m%d')

    # Construct log file name
    log_file = f"{date_str}.{logfilename}"
    log_file += ".log"

    # Ensure the log file is in the specified output directory
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists
        log_file = os.path.join(output_dir, log_file)

    # Initialize the logger
    logger = logging.getLogger(logfilename)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # Create file handler to log to a file
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)

    # Create console handler to print to console
    ch = logging.StreamHandler()
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # Formatter for the log messages
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # Add both handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

# --- Function to run shell commands ---
def run_cmd(cmd, verbose=False, logger=None):
    """
    Run a shell command and handle errors.
    Arguments:
        cmd (str): The shell command to run.
        verbose (bool): If True, print the command before running it.
    Returns:
        None
    """

    # Only log the command if logger is set to DEBUG level
    if logger and logger.isEnabledFor(logging.DEBUG):
        logger.debug(f"Running command: {cmd}")
    
    # Run the command using subprocess
    result = subprocess.run(cmd, shell=True)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}", file=sys.stderr)
        sys.exit(1)

# --- Function to parse command-line arguments ---
def parse_args():
    """
    Set up the argument parser for the script.

    Arguments:
        None
    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """

    parser = argparse.ArgumentParser(description=f'''
    + {VERSION_NAME} v{VERSION} +
    
    Example usage:
        python cellsnp_vcf_merger.py --input-dir /path/to/input --output-dir /path/to/output --study-name study_name
    ''',
    epilog=f'''
    + {VERSION_NAME} v{VERSION}. {COPYRIGHT} \n{COPYRIGHT_TEXT}+''', 
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--input-dir', '-i', type=str, required=True, help='Path to the input directory. Required.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--study-name', '-s', type=str, help='Study name (matches directory within --input-dir). Required unless --merge-all-studies is used.')
    group.add_argument('--merge-all-studies', '-m', action='store_true', help='If set, will merge all studies in the input directory instead of a single study.')
    parser.add_argument('--output-dir', '-o', type=str, required=False, help='Path to the output directory. Optional.')
    parser.add_argument('--clean-up', '-c', action='store_true', help='Remove intermediate merged VCF. Optional.')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print extra information. Optional.')
    parser.add_argument('--debug', '-d', action='store_true', help='Print debug information. Note: this creates a lot of output - think carefully. Optional.')
    parser.add_argument('--version', '-V', action='version', version=f'%(prog)s {VERSION} ({VERSION_DATE}).')
    return parser.parse_args()

# --- Function to merge all studies ---
def merge_all_studies(args):
    base_dir = Path(args.input_dir).resolve()
    if not base_dir.exists():
        logger.error(f"Input directory does not exist: {base_dir}")
        sys.exit(1)

    # Setting up the output directory for merged VCFs
    merged_output_dir = Path(args.output_dir or base_dir / "MetaPlaq2_cellSNP_AllStudies").resolve()
    merged_output_dir.mkdir(parents=True, exist_ok=True)

    # Collect all VCF files from all studies
    vcf_files = []
    for study_path in sorted(base_dir.glob("*")):
        if study_path.is_dir():
            study_name = study_path.name
            merged_vcf = study_path / f"{study_name}.merged.sorted.vcf.gz"
            if merged_vcf.exists():
                vcf_files.append(str(merged_vcf))
            else:
                logger.warning(f"Merged VCF not found: {merged_vcf}")

    if not vcf_files:
        logger.error(f"No per-study merged VCF files found in {base_dir}. Cannot merge.")
        sys.exit(1)
    
    logger.info(f"Found {len(vcf_files)} per-study merged VCF files across all studies.")

    merged_vcf = merged_output_dir / "MetaPlaqv2.cellSNP.AllStudies.merged.vcf"
    sorted_vcf = merged_output_dir / "MetaPlaqv2.cellSNP.AllStudies.merged.sorted.vcf.gz"

    logger.info("Merging all per-study merged VCF files into a single VCF...")
    input_list = " ".join(f'"{str(f)}"' for f in vcf_files)
    with open(merged_output_dir / "AllStudies.input_vcfs.txt", "w") as f:
        f.write("\n".join(map(str, vcf_files)))

    run_cmd(f"bcftools merge -O v -o {merged_vcf} {input_list}", verbose=args.verbose, logger=logger)
    logger.info("Sorting the merged VCF...")
    run_cmd(f"bcftools sort -O z -o {sorted_vcf} {merged_vcf}", verbose=args.verbose, logger=logger)
    logger.info("Indexing the sorted VCF file...")
    run_cmd(f"bcftools index -t {sorted_vcf}", verbose=args.verbose, logger=logger)

    if args.clean_up:
        if merged_vcf.exists():
            if args.verbose:
                logger.info(f"Removing intermediate merged VCF: {merged_vcf}")
            merged_vcf.unlink(missing_ok=True)

    logger.info(f"[DONE] Final sorted VCF file for all studies: {sorted_vcf}")

# --- Main function ---
def main():
    """
    Main function to merge CellSNP VCFs by study.

    This script merges all *_cellSNP_sorted.cells.vcf.gz files from a given study directory,
    then sorts and indexes the merged VCF using bcftools.
    It requires the bcftools command-line tool to be installed and available in the PATH.

    It can also merge all studies in the input directory if the --merge-all-studies flag is set.
    It logs the process and can optionally clean up intermediate files.
    It also provides command-line arguments for customization, including input and output directories,
    study name, verbosity, and debug mode.
    
    Usage:
        python cellsnp_vcf_merger.py --input-dir /path/to/input --output-dir /path/to/output --study-name study_name [--verbose] [--debug]
    
    Arguments:
        --input-dir, -i: Path to the input directory containing the study directories. Required.
        --output-dir, -o: Path to the output directory where the merged VCF will be saved. Optional.
        --study-name, -s: Name of the study (matches directory within --input-dir). Required.
        --clean-up, -c: Remove intermediate merged VCF after sorting and indexing. Optional.
        --merge-all-studies, -m: If set, will merge all studies in the input directory instead of a single study. Optional.
        --verbose, -v: Print extra information. Optional.
        --debug, -d: Print debug information. Note: this creates a lot of output - think carefully. Optional.
        --version, -V: Show version information and exit.
    Returns:
        None
    """
       
    # Parse command-line arguments
    args = parse_args()
    
    # Start the timer
    start_time = time.time()
    
    # Get today's date for logging
    today_date = datetime.now()
    today_str = today_date.strftime("%Y%m%d")

    # Get the project directory -- this is one level down from the script's location
    # This assumes the script is located in the 'cellsnp' directory
    project_dir = Path(__file__).resolve().parent.parent
        # Setup logger after parsing arguments
    global logger
    # Ensure the log is written to the cellsnp directory
    cellsnp = project_dir / 'cellsnp'
    logger = setup_logger(logfilename='cellsnp_vcf_merger', verbose=args.verbose, debug=args.debug, output_dir=cellsnp)

    # Start the script
    logger.info(f"+ {VERSION_NAME} v{VERSION} ({VERSION_DATE}) +")
    logger.info(f"\nStarting extraction job {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")

    # Set up directories
    # Resolve the absolute base input directory (where all studies live)
    base_dir = Path(args.input_dir).resolve()

    # Only define study_dir if not merging all studies; ensure study_name is only accessed when defined
    if args.merge_all_studies:
        study_name_used = "AllStudies"
        study_dir = Path(args.output_dir).resolve() if args.output_dir else (base_dir / ".AllStudies.")
    else:
        if not args.study_name:
            logger.error("Either --study-name must be provided or --merge-all-studies must be enabled.")
            sys.exit(1)
        study_name_used = args.study_name
        study_dir = base_dir / args.study_name
        if not study_dir.exists():
            logger.error(f"Study directory does not exist: {study_dir}. Please check the input directory and study name.")
            print(f"[ERROR] Study directory does not exist: {study_dir}", file=sys.stderr)
            sys.exit(1)
    
    # Set the output directory (either user-specified or default to study_dir)
    output_dir = Path(args.output_dir).resolve() if args.output_dir else study_dir
    if args.output_dir and not str(output_dir).startswith(str(study_dir)):
        logger.warning(f"Output directory is outside the study directory: {output_dir}")

    # Ensure the output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Print extra information 
    logger.info(f"Input directory.........: {base_dir}")
    if args.output_dir:
        logger.info(f"Output directory........: {Path(args.output_dir).resolve()}")
    else:
        logger.info(f"Output directory........: {study_dir}")
    logger.info(f"Study directory.........: {study_dir}")
    logger.info(f"Study name..............: {study_name_used}")
    logger.info(f"Clean up intermediate...: {'Yes' if args.clean_up else 'No (default)'}")
    logger.info(f"Merge all studies.......: {'Yes' if args.merge_all_studies else 'No (default)'}")
    logger.info(f"Verbose mode............: {'On' if args.verbose else 'Off (default)'}")
    logger.info(f"Debug mode..............: {'On' if args.debug else 'Off (default)'}")
    logger.info(f"Running version.........: v{VERSION} ({VERSION_DATE})")
    logger.debug(f"Project directory.......: {project_dir}")
    logger.debug(f"Cellsnp directory.......: {cellsnp}")
    logger.debug(f"Script location.........: {Path(__file__).resolve()}")
    logger.info(f"Python version..........: {sys.version}")
    logger.info(f"Current date............: {today_date}")
    logger.info(f"Current time............: {datetime.now().strftime('%H:%M:%S')}\n")

    # Check if bcftools is installed
    import shutil
    if shutil.which("bcftools") is None:
        logger.error("bcftools not found in PATH. Please install it before running this script.")
        print("[ERROR] bcftools not found in PATH. Please install it before running this script.", file=sys.stderr)
        sys.exit(1)

    # If merging all studies is requested, call the merge_all_studies function
    if args.merge_all_studies:
        logger.info("Merging all studies in the input directory...")
        merge_all_studies(args)
        logger.info("All studies merged successfully.")
        return

    # If merging a single study, proceed with the main logic
    logger.info("Merging VCF files for a single study...")
    # Get the list of VCF files in the study directory
    vcf_files = sorted(study_dir.rglob("*_cellSNP_sorted.cells.vcf.gz"))
    if not vcf_files:
        # If no VCF files are found, print an error and exit
        logger.error(f"No VCF files found in {study_dir}. Please check the directory and ensure it contains the expected VCF files.")
        print(f"[ERROR] No VCF files found in {study_dir}", file=sys.stderr)
        sys.exit(1)
    logger.info(f"Found {len(vcf_files)} VCF files in {study_dir}.")
    logger.debug(f"VCF files: {', '.join(str(f) for f in vcf_files)}")
    
    # Create temp directory for reheadered VCFs
    tmp_dir = output_dir / "tmp_renamed_vcfs"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    renamed_vcf_files = []

    logger.info("Renaming VCF samples using bcftools reheader...")

    for vcf in vcf_files:
        sample_name = f"{study_name_used}_{vcf.parent.name.replace('/', '_')}"
        rename_txt = tmp_dir / f"{sample_name}_rename.txt"
        renamed_vcf = tmp_dir / f"{sample_name}.vcf.gz"
        
        # Write renaming map: assume original sample is 'Sample_0'
        with open(rename_txt, "w") as f:
            f.write(f"Sample_0\t{sample_name}\n")
        
        # Run bcftools reheader
        logger.debug(f"Renaming sample in header {vcf} to {renamed_vcf} using {rename_txt}")
        run_cmd(f"bcftools reheader -s {rename_txt} -o {renamed_vcf} {vcf}", verbose=args.verbose, logger=logger)
        # Index the renamed VCF
        run_cmd(f"bcftools index -t {renamed_vcf}", verbose=args.verbose, logger=logger)
        if args.verbose:
            logger.info(f"Verbose mode -- Renamed VCF: {renamed_vcf} (original: {vcf})")
        renamed_vcf_files.append(renamed_vcf)

    vcf_files = renamed_vcf_files
    if args.verbose:
            logger.info(f"Verbose mode -- Using reheadered VCFs from: {tmp_dir}")

    # Create the merged VCF file name and sorted VCF file name
    merged_vcf = output_dir / f"{study_name_used}.merged.vcf"
    sorted_vcf = output_dir / f"{study_name_used}.merged.sorted.vcf.gz"

    if args.verbose:
            logger.info(f"Verbose mode -- Merging VCF files into: {merged_vcf}")
            logger.info(f"Verbose mode -- Sorted VCF will be saved as: {sorted_vcf}")

    # Create the input list for bcftools merge
    input_list = " ".join(f'"{str(f)}"' for f in vcf_files)
    # Save the input list to a file
    with open(output_dir / f"{study_name_used}.input_vcfs.txt", "w") as f:
        f.write("\n".join(map(str, vcf_files)))
    logger.debug(f"Input list for bcftools merge: {input_list}")
    # Run bcftools merge, sort, and index
    logger.info(f"Merging VCF files...")
    run_cmd(f"bcftools merge -O v -o {merged_vcf} {input_list}", verbose=args.verbose, logger=logger)
    logger.info(f"Sorting the merged VCF...")
    run_cmd(f"bcftools sort -O z -o {sorted_vcf} {merged_vcf}", verbose=args.verbose, logger=logger)
    logger.info(f"Indexing the sorted VCF file...")
    run_cmd(f"bcftools index -t {sorted_vcf}", verbose=args.verbose, logger=logger)

    # Optionally clean up. If you want to remove the merged VCF file after sorting and indexing
    if merged_vcf.exists() and args.clean_up:
        logger.info(f"Cleaning up: removing the intermediate merged VCF file: {merged_vcf}")
        # Remove the merged VCF file
        merged_vcf.unlink(missing_ok=True)
    else:
        logger.info(f"Keeping the intermediate merged VCF file: {merged_vcf}")
    # Print the final sorted VCF file
    logger.info(f"[DONE] Final sorted VCF file: {sorted_vcf}")

    # Calculate and print execution time
    elapsed_time = time.time() - start_time
    time_delta = timedelta(seconds=elapsed_time)
    formatted_time = str(time_delta).split('.')[0]
    logger.info(f"Script executed on {today_date.strftime('%Y-%m-%d')}. Total execution time: {formatted_time}.")

    # Clean up temporary directory with renamed VCFs
    if args.clean_up and tmp_dir.exists():
        import shutil
        logger.info(f"Cleaning up: removing temporary directory with reheadered VCFs: {tmp_dir}")
        shutil.rmtree(tmp_dir)
    else:
        logger.info(f"Keeping temporary directory with reheadered VCFs: {tmp_dir}")
    
    # Print the version and license information
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}). {COPYRIGHT}")
    # logger.info(COPYRIGHT_TEXT.strip())

if __name__ == "__main__":
    main()
# End of the script