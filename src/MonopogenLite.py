#!/usr/bin/env python3
"""
The main interface of MonopgenLite.

MonopogenLite is a lightweight version of Monopogen, a tool for single-cell variant calling and phasing.
It is designed to perform germline variant calling and phasing from single-cell sequencing data.
It is a simplified version of Monopogen, focusing on germline variant calling and phasing, since 
de novo calling of variants is not feasible with most single-cell sequencing data. 

Arguments:
    - `--app_path`: Path to the application directory containing the Beagle JAR file.
    - `--nthreads`: Number of threads to use for processing.
    - `--debug`: Enable debug mode for verbose output.
    - `--verbose`: Enable verbose output.
    - `--out`: Output directory where results will be saved.
    - `--region`: Path to a file containing regions to process, one per line.
    - `--reference`: Path to the reference genome FASTA file.
    - `--imputation_panel`: Path to the imputation panel VCF file.
    - `--max_softClipped`: Maximum number of soft-clipped bases allowed in reads.
    - `--max_mismatch`: Maximum number of mismatches allowed in reads (for pre-processing).
    - `--platform_library`: Platform library used for sequencing (e.g., Illumina, PacBio).
    - `--min_read_length`: Minimum read length for filtering reads (for pre-processing).
    - `--umi_collapse`: UMI collapse option, can be a boolean or a string indicating the UMI tag.
Returns:
    - pre-processed BAM files in the specified output directory needed for germline variant calling.
    - Germline variant calling and phasing results in the specified output directory.
    - Job scripts for each region processed, which can be run in parallel.

This script is part of the MonopogenLite project, which is licensed under the MIT License.

"""

# Importing the required libraries -- 2024-09-17
# import general libraries
import argparse # parser for command-line options, arguments and sub-commands
import sys # system-specific parameters and functions
import os # operating system dependent functionality
import logging # logging module
import shutil # high-level file operations
import glob # Unix style pathname pattern expansion
import re # regular expressions
import time # handle time
from datetime import datetime
import subprocess # run shell commands
import multiprocessing as mp
from multiprocessing import Pool
from pathlib import Path
import pysam

# import data manipulation and analysis libraries
import pandas as pd # data manipulation and analysis
import numpy as np # numerical computing
import gzip # read and write gzip files

# import libraries to handle BAM and VCF files
import pysam # handle BAM files
from pysam import VariantFile # handle VCF files

# import submodules
# from bamProcess import * 
from germline import *

# Version and license information
VERSION_NAME = 'MonopogenLite'
VERSION = '1.4.7'
VERSION_DATE = '2026-11-09'
COPYRIGHT = 'Copyright 1979-2026. Jinzhuang Dou | jdou1 [at] mdanderson [dot] org; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
COPYRIGHT_TEXT = '''
The MIT License (MIT).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
'''

# Setting the library path -- 2024-08-08
LIB_PATH = Path(__file__).resolve().parent / "pipelines/lib"
LIB_PATH = str(LIB_PATH)
if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)

PIPELINE_BASEDIR = str(Path(sys.argv[0]).resolve().parent)
CFG_DIR = str(Path(PIPELINE_BASEDIR) / "cfg")

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)

# Function to optionally add file logging
def setup_logger(logfile_path, level=logging.DEBUG):
    """
    Set up the logger to log to a file and console.
    Args:
        logfile_path (str): Path to the log file.
        level (int): Logging level (default: logging.DEBUG).
    """

    handler = logging.FileHandler(logfile_path)
    handler.setFormatter(logging.Formatter('[{asctime}] {levelname:8s} {filename} {message}', style='{'))
    logger.addHandler(handler)
    logger.setLevel(level)

# Function to print errors if any exist -- 2024-08-08
def error_check(all, output, step):
    """
    Check if all expected jobs have been completed successfully.
    Args:
        all (list): List of all expected job IDs.
        output (list): List of completed job IDs.
        step (str): The step of the pipeline being checked.
    """

    job_fail = 0
    for id in all:
        if id not in output:
            logger.error("In "+ step + " step " + id + " failed!")
            job_fail = job_fail + 1

    if job_fail > 0:
        logger.error("Failed! See instructions above.")
        sys.exit(1)

# Utility function to log tool versions -- 2024-06-05
def log_tool_versions(logger):
    """
    Log the versions of essential tools used in the pipeline.
    Args:
        logger (logging.Logger): Logger instance to log the tool versions.
    """

    # Use the global 'args' variable if present, else skip beagle version
    try:
        app_path = args.app_path
    except Exception:
        app_path = None
    tools = {
        "bcftools": ["bcftools", "--version"],
    }
    if app_path:
        tools["beagle"] = ["java", "-jar", str(Path(app_path) / "beagle.jar")]
    for name, cmd in tools.items():
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            version_info = result.stdout.strip().split('\n')[0] or result.stderr.strip().split('\n')[0]
            logger.info(f"{name} version: {version_info}")
        except Exception as e:
            logger.warning(f"Failed to retrieve version for {name}: {e}")

# Function to build sample commands for variant calling and phasing
def build_sample_commands(samples, chrom, out_path, args):
    commands = []
    for sample in samples:
        bam_path = sample["bam"]
        sampleID = sample["sampleID"]
        output_vcf = out_path / "VCF" / f"{sampleID}.chr{chrom}.bcf"
        # Updated: pass bam_path, chrom, args.reference, out_path
        bcftools_cmd = generate_bcftools_command(bam_path, chrom, args.reference, out_path)
        # bcftools_cmd = generate_bcftools_command(bam_path, sampleID, chrom, out_path)
        if hasattr(args, 'impute') and args.impute:
            if hasattr(args, 'genotype') and args.genotype:
                beagle_cmd = generate_beagle_cmd_gt(output_vcf, args.app_path, out_path)
            else:
                beagle_cmd = generate_beagle_cmd_gp(output_vcf, args.app_path, out_path)
            command = f"{bcftools_cmd} && {beagle_cmd}"
        else:
            command = bcftools_cmd
        commands.append((sampleID, command))
    return commands

# Function to write a job script for germline variant calling
def write_job_script(jobid, command, out_path, version=VERSION, verbose=False):
    script_path = out_path / f"runGermline_{jobid}.sh"
    vcf_file_name = f"{out_path}/monopogen_{jobid}"
    vcf_extension = ".vcf.gz"
    # pipefail and '|| exit 1' make the script exit with an error when any step fails, so the failure is
    # reported instead of the script ending successfully; sorting/indexing is part of the command -- 2026-11-09
    slurm_script_content = f"""#!/bin/bash
echo "[MonopogenLite germline.py v{version}] Start time: $(date)"
echo "[MonopogenLite germline.py v{version}] Running job: {jobid}"

source ~/.bashrc
micromamba activate monopogen
set -o pipefail

{command} || {{ echo "[MonopogenLite germline.py v{version}] Job {jobid} FAILED: $(date)"; exit 1; }}

#rm {vcf_file_name}{vcf_extension} {vcf_file_name}.gp{vcf_extension} {vcf_file_name}.phased{vcf_extension}

echo "[MonopogenLite germline.py v{version}] End time: $(date)"
"""
    write_slurm_script(script_path, slurm_script_content, verbose=verbose)
    return str(script_path)

# Helper function to write SLURM (shell) scripts
def write_slurm_script(path, content, verbose=False):
    """
    Write a SLURM script to the specified path.
    Args:
        path (str or Path): The path where the script will be written.
        content (str): The content of the SLURM script.
        verbose (bool): If True, print a message indicating the script was written.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(content)
    # Make the script executable.
    os.system(f'chmod +x {path}')
    if verbose:
        print(f"  - SLURM script written to: {path}")

# Function to generate bcftools command for mpileup and filtering
def generate_bcftools_command(bam_filter, jobid, reference, out):
    """
    Generate the bcftools command for mpileup and filtering.
    Args:
        bam_filter (str): Path to the list of BAM files.
        jobid (str): Job ID (which is the chrom[osome] number) or region to process.
        reference (str): Path to the reference genome FASTA file.
        out (str): Output directory where the VCF file will be saved.
    """
    # FIXME: mpileup only computes the PL. Added bcftools call this actually computes the genotypes based on the PL.
    return (
        f"{bcftools} mpileup -b {bam_filter} --fasta-ref {reference} --regions {jobid} --min-MQ 20 --min-BQ 20 --annotate FORMAT/DP "
        f"| {bcftools} call -mv -Ou  " 
        f"| {bcftools} norm -m-both --rm-dup both --check-ref wx --fasta-ref {reference} "
        f"| {bcftools} annotate --set-id '%CHROM:%POS:%REF:%ALT' "
        f"| {bcftools} filter --exclude 'ALT !~ \"^[ATGC]$\"' "
        f"| {bcftools} +fill-tags -Oz -o {out}/monopogen_{jobid}.vcf.gz")

# FIXME: 24-01-2026: Adjusted paths of gt and gp to ensure phasing happens correctly! (Damian)
# Function to generate Beagle command for genotype probabilities
def generate_beagle_cmd_gp(jobid, out, imputation_vcf, chrom, nthreads_downsample):
    """
    Generate the Beagle command for genotype probabilities.
    Args:
        jobid (str): Job ID or region cd to process.
        out (str): Output directory where the genotype probabilities will be saved.
        imputation_vcf (str): Path to the imputation panel VCF file.
        chrom (str): Chromosome to process.
        nthreads_downsample (int): Number of threads to use for downsampling.
    """

    return (
        f"{java} -Xmx20g -jar {beagle} gl={out}/monopogen_{jobid}.vcf.gz ref={imputation_vcf} "
        f"chrom={chrom} out={out}/monopogen_{jobid}.gp impute=false modelscale=2 "
        f"nthreads={nthreads_downsample} gprobs=true niterations=0"
    )

# Function to generate Beagle command for genotype phasing
# FIXME: 24-01-2026: Adjusted paths of gt and gp to ensure phasing happens correctly! (Damian)
def generate_beagle_cmd_gt(jobid, out, imputation_vcf, chrom, nthreads_downsample):
    """
    Generate the Beagle command for genotype phasing.
    Args:
        jobid (str): Job ID or region to process.
        out (str): Output directory where the phased genotypes will be saved.
        imputation_vcf (str): Path to the imputation panel VCF file.
        chrom (str): Chromosome to process.
        nthreads_downsample (int): Number of threads to use for downsampling.
    """

    return (
        f"{java} -Xmx20g -jar {beagle} gt={out}/monopogen_{jobid}.gp.vcf.gz ref={imputation_vcf} "
        f"chrom={chrom} out={out}/monopogen_{jobid}.phased impute=false modelscale=2 "
        f"nthreads={nthreads_downsample} gprobs=true niterations=0"
    )

# Function to make output directories
def prepare_output_dirs(base_out, subdirs, verbose=False):
    """
    Create output directories if they do not exist.
    Args:
        base_out (str or Path): Base output directory where subdirectories will be created.
        subdirs (list): List of subdirectory names to create.
        verbose (bool): If True, print messages about created directories.
    """
    base_out = Path(base_out)
    base_out.mkdir(parents=True, exist_ok=True)
    for subdir in subdirs:
        path = base_out / subdir
        path.mkdir(parents=True, exist_ok=True)
        if verbose:
            print(f"  - Created directory: {path}")

# Function to read a sample list based on the region file and return a list of dicts with a "bam" key
def read_sample_list_file(region_file_path, out):
    """
    Read the sample list based on the region file.

    Args:
        region_file_path (str): Path to the region file.
        out (str): Output directory.

    Returns:
        list of dict: A list of samples where each sample is a dictionary with at least a 'bam' key.
    """
    sample_set = set()
    with open(region_file_path) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) >= 1:
                chrom = parts[0]
                # FIXME: Adjust the file path to the correct folder with the lists containing the BAMS per chromosome.
                bam_file_path = Path(out) / ".." / "ListedBamsPerChr" / f"{chrom}.lst"
                print(
                    f"  -- DEBUGGING: The BAM file path is: [{bam_file_path}].")
                with bam_file_path.open() as bf:
                    for line in bf:
                        bam = line.strip()
                        sample_set.add(bam)

    samples = [{"bam": bam} for bam in sorted(sample_set)]
    return samples

# Chromosomes MonopogenLite can call; the imputation panel covers 1-22 and X -- 2026-11-09
SUPPORTED_CHROMOSOMES = [str(n) for n in range(1, 23)] + ["X"]

def normalize_chromosome(name):
    """
    Return the chromosome as '1'-'22' or 'X' (without 'chr'), or None if it is not supported.
    Args:
        name (str or int): Chromosome, e.g. 1, '1', 'chr1', 'X', 'chrX'.
    """
    bare = str(name).strip()
    if bare.lower().startswith("chr"):
        bare = bare[3:]
    bare = bare.upper()
    return bare if bare in SUPPORTED_CHROMOSOMES else None

def read_regions_for_chromosome(region_file, chromosome):
    """
    Return the regions of one chromosome in the region file, as 'chrN' or 'chrN:start-end'.
    Lines are 'chromosome' or 'chromosome,start,end'; blank lines and lines starting with '#' are skipped.
    Args:
        region_file (str): Path to the region file.
        chromosome (str): Chromosome as returned by normalize_chromosome().
    """
    regions = []
    with open(region_file) as f_in:
        for line_number, line in enumerate(f_in, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            record = [field.strip() for field in line.split(",")]
            if len(record) not in (1, 3):
                print(f"ERROR: Line {line_number} of the region file [{region_file}] should have 1 column (chromosome) "
                      f"or 3 comma-delimited columns (chromosome,start,end): [{line}].")
                sys.exit(1)
            if normalize_chromosome(record[0]) != chromosome:
                continue
            if len(record) == 1:
                regions.append(f"chr{chromosome}")
            else:
                regions.append(f"chr{chromosome}:{record[1]}-{record[2]}")
    return regions

def build_chromosome_bam_list(out, preprocess_dir, chrom_label):
    """
    Return the path to OUT/ListedBamsPerChr/<chrom_label>.lst and the number of BAMs in it.
    If the list does not exist, it is built from the preProcess output: every
    <preprocess_dir>/*/Bam/*_<chrom_label>.filter.bam. Every listed BAM must exist and be indexed.
    Args:
        out (str or Path): Output directory of germline.
        preprocess_dir (str or Path): Directory holding one preProcess output folder per sample.
        chrom_label (str): Chromosome with 'chr' prefix, e.g. 'chr1'.
    """
    list_path = Path(out) / "ListedBamsPerChr" / f"{chrom_label}.lst"
    found = sorted(str(p) for p in Path(preprocess_dir).glob(f"*/Bam/*_{chrom_label}.filter.bam"))
    if list_path.is_file():
        listed = [line.strip() for line in list_path.open() if line.strip()]
        logger.info(f"Using the existing BAM list [{list_path}] with {len(listed)} BAM(s).")
        if sorted(listed) != found:
            logger.warning(f"The existing BAM list [{list_path}] differs from the {len(found)} filtered BAM(s) found in "
                           f"[{preprocess_dir}/*/Bam/]. Delete the list to let germline rebuild it.")
    else:
        if not found:
            print(f"ERROR: No filtered BAM files (*_{chrom_label}.filter.bam) found in [{preprocess_dir}/*/Bam/]. "
                  f"Run preProcess first, or set --preprocess-dir.")
            sys.exit(1)
        list_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = list_path.with_name(list_path.name + ".tmp")
        tmp_path.write_text("\n".join(found) + "\n")
        os.replace(tmp_path, list_path)
        listed = found
        logger.info(f"Wrote the BAM list [{list_path}] with {len(found)} BAM(s) found in [{preprocess_dir}/*/Bam/].")
    missing = [bam for bam in listed if not Path(bam).is_file() or not Path(bam + ".bai").is_file()]
    if missing:
        print(f"ERROR: {len(missing)} BAM(s) in [{list_path}] are missing or not indexed, e.g. {missing[:3]}.")
        sys.exit(1)
    return list_path, len(listed)

# Function to perform germline variant calling -- 2024-08-08
def germline(args):
    """
    Perform germline variant calling and phasing from single-cell sequencing data.
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """

    logger.info("Performing germline variant calling with the following arguments.")
    print(f"  -- DEBUGGING: Germline before print_parameters_given()")
    print_parameters_given(args)

    logger.info("Checking existence of essenstial resource files.")
    print(f"  -- DEBUGGING: Germline before validate_user_setting()")
    validate_user_setting_germline(args)

    logger.info("Checking dependencies.")
    print(f"  -- DEBUGGING: Germline before check_dependencies()")
    check_dependencies(args)

    # Create necessary directories -- 2024-09-17
    if not args.out:
        print(f"ERROR: Output directory not specified!")
        sys.exit(1)
    # No longer neccessary for per region variant calling.
    #prepare_output_dirs(args.out, ['germline', 'scripts'], verbose=args.verbose)
    out_path = Path(args.out)

    # Load samples from the region file (NOT RELEVANT FOR REGION SPECIFIC VARIANT CALLING)
    #samples = read_sample_list_file(region_file_path=args.region, out=args.out)
    # Convert to expected sample dict format
    #samples = [{"sampleID": Path(bam["bam"]).stem.replace(".filter", "").replace(".bam", ""), "bam": bam["bam"]} for bam in samples]
    #if args.debug:
        #print(f"  -- DEBUGGING: Sample list loaded with {len(samples)} samples.")
        #print(f"  -- DEBUGGING: Sample list: {samples}")

    # Select the regions of this chromosome from the region file -- 2026-11-09
    # Regions are matched on chromosome name instead of line number, so the region file may list the
    # chromosomes in any order, with or without the 'chr' prefix, and may hold several start-end regions
    # per chromosome (one job each).
    chromosome = normalize_chromosome(args.chromosome)
    if chromosome is None:
        print(f"ERROR: --chromosome '{args.chromosome}' is not supported; use 1-22 or X, with or without 'chr'.")
        sys.exit(1)
    chrom_label = f"chr{chromosome}"
    if args.verbose:
        print(f"  - Checking the region file [{args.region}] for {chrom_label}.")
    regions = read_regions_for_chromosome(args.region, chromosome)
    if not regions:
        print(f"ERROR: The region file [{args.region}] has no line for {chrom_label}.")
        sys.exit(1)
    logger.info(f"Found {len(regions)} region(s) for {chrom_label}: {', '.join(regions)}.")

    # List of filtered BAMs (one per sample) for this chromosome, built from the preProcess output -- 2026-11-09
    preprocess_dir = args.preprocess_dir or args.out
    bam_filter, N_sample = build_chromosome_bam_list(out_path, preprocess_dir, chrom_label)
    if args.verbose:
        print(f"  - {N_sample} sample BAM(s) listed in [{bam_filter}].")

    # NEW code -- 2024-09-17
    # check the imputation panel file
    if args.verbose:
        print(f"  - Checking the imputation panel file: [{args.imputation_panel}].")
    imputation_vcf = Path(args.imputation_panel) / f"1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.{chrom_label}.vcf.gz"
    if not imputation_vcf.is_file():
        print(f"ERROR: The imputation panel file [{imputation_vcf}] cannot be found!")
        sys.exit(1)

    # Beagle gets half of the threads, but at least one -- 2026-11-09
    nthreads_downsample = max(1, args.nthreads // 2)
    if args.debug:
        print(f"  -- DEBUGGING: number of threads used: {args.nthreads}; Beagle uses {nthreads_downsample}.")

    joblst = []
    for jobid in regions:
            if args.verbose:
                print(f"  - Starting germline variant calling for region {jobid}.")

            # NEW CODE with bcftools -- 2025-06-05
            # germline variant calling from mpileup
            # https://samtools.github.io/bcftools/bcftools.html
            # https://www.biostars.org/p/425139/
            # https://www.biostars.org/p/418738/
            # mpileup a single region
            # -b list of input BAM files
            # --regions region to include
            # --min-MQ Minimum mapping quality for an alignment to be used [0]
            # --min-BQ Minimum base quality for a base to be considered [13]
            # --annotate FORMAT/DP
            # norm -m-both --rm-dup both > split multi-allelics, both INDEL and SNP, into bi-allelic, but remove duplicate alleles based on chrom, pos, ref, alt
            # --check-ref wx -f > check reference allele (-w warn only; -x fix issues) and use reference data (-f)
            # what to do when incorrect or missing REF allele is encountered:
            # exit (e), warn (w), exclude (x), or set/fix (s) bad sites.
            # The w option can be combined with x and s.
            # Note that s can swap alleles and will update genotypes (GT) and AC counts,
            # but will not attempt to fix PL or other fields.
            # Also note, and this cannot be stressed enough,
            # that s will NOT fix strand issues in your VCF, do NOT use it for that purpose!!!
            # (Instead see http://samtools.github.io/bcftools/howtos/plugin.af-dist.html and <https://samtools.github.io/bcftools/howtos/plugin.fixref.html>.)
            # annotate --set-id '%CHROM:%POS:%REF:%ALT' > set ID to CHROM:POS:REF:ALT
            # filter -e \'ALT !~ "^[ATGC]$"\ > to remove non-ATGC bases
            # +fill-tags -Oz > fill tag (the INFO field) and compress output

            # the full command
            # bcftools mpileup -b /path/to/MonopogenLite/example/Bam/chr22.filter.bam.lst --check-ref wx --fasta-ref /Users/slaan3/PLINK/references/fasta/hg38.fa --regions chr22 --min-MQ 20 --min-BQ 20 --annotate FORMAT/DP
            # | bcftools norm -m-both --rm-dup both --check-ref wx --fasta-ref /path/to/references/fasta/hg38.fa
            # | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'
            # | bcftools filter --exclude 'ALT !~ "^[ATCG]$"'
            # | bcftools +fill-tags -Oz -o /path/to/MonopogenLite/example/germline_test/chr22.gl.vcf.gz
            cmd1 = generate_bcftools_command(bam_filter, jobid, args.reference, str(out_path))

            # NEW code -- 2024-09-17
            if args.verbose:
                print(f"    * germline variant calling")
            if args.debug:
                print(f"    -- DEBUGGING: Command to run: {cmd1}")

            # NEW code -- 2025-06-05
            # germline variant phasing of genotype probabilities
            cmd3 = generate_beagle_cmd_gp(jobid, str(out_path), imputation_vcf, chrom_label, nthreads_downsample)

            # NEW code -- 2025-06-05
            if args.verbose:
                print(f"    * germline variant phasing of genotype probabilities")
            if args.debug:
                print(f"    -- DEBUGGING: Command to run: {cmd3}")

            # NEW code -- 2025-06-05
            # germline variant phasing of genotypes
            # cmd5 = java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  jobid + ".germline.vcf"  +  " ref=" +  imputation_vcf  +  "  chrom=" + record[0]  + " out="   +  out + "/germline/" + jobid+ ".phased " + "impute=false modelscale=2 nthreads=" + nthreads_downsample + " gprobs=true niterations=0"
            # cmd5 = cmd5 + "\n" + "rm -v " +  out + "/germline/" +  jobid + ".germline.vcf"
            cmd5 = generate_beagle_cmd_gt(jobid, str(out_path), imputation_vcf, chrom_label, nthreads_downsample)

            # NEW code -- 2024-09-17
            if args.verbose:
                print(f"    * germline variant phasing of genotypes")
            if args.debug:
                print(f"    -- DEBUGGING: Command to run: {cmd5}")

            # NEW code -- 2025-06-05
            # Use build_sample_commands and write_job_script to generate job scripts
            # These should replace the above job script generation logic.
            # Assumes build_sample_commands(samples, chr, out_path, args) returns a list of (jobid, command) tuples.
            # Since the current function is per-region, we need to adapt the context.
            # Here, we use the current jobid and command composition as a single job, so we simulate the API.
            # If build_sample_commands and write_job_script are available, use them here.
            # (If not, they need to be integrated from external code.)
            # OBSOLETE code -- 2024-09-17
            # The genotype field typically indicates whether a variant is present in an individual. For example:
            # - 0/0 means the individual is homozygous for the reference allele (no variation).
            # - 0/1 means the individual is heterozygous, carrying one reference allele and one alternate allele.
            # - 1/1 means the individual is homozygous for the alternate allele.
            # Filtering out 0/0 calls means you are removing positions where the individual has
            # the reference allele on both chromosomes and no variant.
            # To increase computational efficiency, we could filter out 0/0 calls from the VCF file
            # when it concerns one sample. However, it is not useful in all cases. For instance,
            # when you are working with multiple samples, you might want to keep the 0/0 calls
            # to identify the reference allele frequency and to perform population genetics analyses.
            # The following command filters out 0/0 calls from the VCF file when it concerns one sample
            # and writes the output to a new VCF file conditional on a flag -z, --homozygote-filter.
            # However, after careful consideration this has changed to always keep the 0/0 calls.
            # If not, we artificially create missing data, where there is none. It is better
            # to filter 0/0 calls later prior to downstream genetic analyses.
            # if N_sample == 1:
            # 	if args.homozygote_filter:
            # 		cmd4 = "zless -S " +  out + "/germline/" + jobid + ".gp.vcf.gz | grep -v  0/0  > " +  out + "/germline/" + jobid + ".germline.vcf"
            # elif N_sample > 1:
            # 	cmd4 = "zless -S " +  out + "/germline/" + jobid + ".gp.vcf.gz   > " +  out + "/germline/" + jobid + ".germline.vcf"
            # f_out.write(cmd4 + "\n")

            # Fixme: Uncommented, because this causes per sample-jobs.
            # commands = build_sample_commands(samples, record[0], out_path, args)
            # for jobid, (_, command) in enumerate(commands):
            # 	write_job_script(jobid, command, out_path, version=VERSION, verbose=args.verbose)
            # for jobid, (_, command) in enumerate(commands):
            # 	script_path = write_job_script(jobid, command, out_path, version=VERSION, verbose=args.verbose)
            # 	joblst.append(f'bash "{script_path}"')


            # Compose the steps requested with --step; one job per region -- 2026-11-09
            # The unphased VCF is sorted and indexed right after calling, so it is kept even if phasing fails.
            unphased_vcf = f"{out_path}/monopogen_{jobid}.vcf.gz"
            sorted_vcf = f"{out_path}/monopogen_{jobid}_sorted.vcf.gz"
            cmd_sort = f"{bcftools} sort {unphased_vcf} -Oz -o {sorted_vcf} && {bcftools} index -t {sorted_vcf}"

            # Beagle stops with an error on a VCF without variants; skip phasing for such a region instead of failing the job.
            def unless_no_variants(cmd):
                return (f'if [ -z "$({bcftools} view -H {unphased_vcf} | head -n 1)" ]; '
                        f'then echo "No variants called in {jobid}; skipping Beagle."; else {cmd}; fi')

            step_commands = {
                "varScan": [cmd1, cmd_sort],
                "varProb": [unless_no_variants(cmd3)],
                "varPhasing": [cmd5],
                "all": [cmd1, cmd_sort, unless_no_variants(f"{cmd3} && {cmd5}")],
            }
            required_input = {"varProb": unphased_vcf,
                              "varPhasing": f"{out_path}/monopogen_{jobid}.gp.vcf.gz"}.get(args.step)
            if required_input and not args.norun and not Path(required_input).is_file():
                print(f"ERROR: --step {args.step} needs [{required_input}], which does not exist. Run the previous step first.")
                sys.exit(1)
            full_cmd = " && ".join(step_commands[args.step])
            if args.debug:
                print(f"  -- DEBUGGING: The full command for region {jobid} (step {args.step}) is\n{full_cmd}.")

            script_path = write_job_script(jobid, full_cmd, out_path, version=VERSION, verbose=args.verbose)
            joblst.append(f'bash "{script_path}"')

    # pool the shell scripts in the job list when the norun flag is not set
    # Each job runs Beagle with nthreads_downsample threads, so run only as many jobs at once as the threads allow.
    if not args.norun:
        with Pool(processes=max(1, min(len(joblst), args.nthreads // nthreads_downsample))) as pool:
            print(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}) \nMonopogenLite: SNV calling and phasing from single-cell sequencing data.")
            print("\nSettings:")
            print(f"  - App path...............: {args.app_path}")
            print(f"  - Threads................: {args.nthreads} (downsample for BEAGLE phasing: {nthreads_downsample})")
            print(f"  - Debug mode.............: {args.debug}")
            print(f"  - Verbosity..............: {args.verbose}")
            print(f"  - Output.................: {args.out}")
            print(f"  - Regions................: {args.region}")
            print(f"  - Reference..............: {args.reference}")
            print(f"  - Imputation panel.......: {args.imputation_panel}")
            print(f"  - Maximum soft-clipped...: {args.max_softClipped}")
            if args.subcommand == "preProcess":
                print(f"  - Maximum mismatch.......: {args.max_mismatch}")
                print(f"  - Platform library.......: {args.platform_library}")
                print(f"  - Minimum read length....: {args.min_read_length}")
                print(f"  - UMI collapse...........: {args.umi_collapse}")
            print(f"\nRunning the germline variant calling pipeline in a pool.")
            print(f"> Jobs submitted:")
            print(joblst)
            result = pool.map(runCMD, joblst)
        # runCMD returns None for a failed job; report it and exit with an error so SLURM marks the job as failed -- 2026-11-09
        failed_jobs = [job for job, res in zip(joblst, result) if res is None]
        if failed_jobs:
            logger.error(f"{len(failed_jobs)} of {len(joblst)} job(s) failed: {failed_jobs}")
            sys.exit(1)
        print(f"\nAll jobs submitted and run.")
        print(f"{VERSION_NAME} v{VERSION}. {COPYRIGHT}\n")
    else:
        print(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}) \nMonopogenLite: SNV calling and phasing from single-cell sequencing data.")
        print("\nSettings:")
        print(f"  - App path...............: {args.app_path}")
        print(f"  - Threads................: {args.nthreads} (downsample for BEAGLE phasing: {nthreads_downsample})")
        print(f"  - Debug mode.............: {args.debug}")
        print(f"  - Verbosity..............: {args.verbose}")
        print(f"  - Output.................: {args.out}")
        print(f"  - Regions................: {args.region}")
        print(f"  - Reference..............: {args.reference}")
        print(f"  - Imputation panel.......: {args.imputation_panel}")
        print(f"  - Maximum soft-clipped...: {args.max_softClipped}")
        if args.subcommand == "preProcess":
            print(f"  - Maximum mismatch.......: {args.max_mismatch}")
            print(f"  - Platform library.......: {args.platform_library}")
            print(f"  - Minimum read length....: {args.min_read_length}")
            print(f"  - UMI collapse...........: {args.umi_collapse}")
        print(f"\nGenerated the job scripts only. The jobs will not be run.")
        print(f"> Jobs generated:")
        print(joblst)
        print(f"\nExiting.")
        print(f"{VERSION_NAME} v{VERSION}. {COPYRIGHT}\n")

# Added by Teresa with help of ClaudeAI (09-11-26).
def resolve_contig(bam, chrom):
    """
    Return the contig name used for a chromosome in a BAM file (e.g. '1' OR 'chr1') or None.
    This function is used in count_and_filter_reads() to resolve how the chromosomes are named in a BAM file.

    Args:
        bam (pysam.AlignmentFile): Open BAM file to be filtered in count_and_filter_reads().
        chrom (str or int): Chromosome label to look up, either with or without 'chr' prefix.
    """
    name = str(chrom)
    bare_name = name[3:] if name.startswith("chr") else name
    refs = set(bam.references)
    for candidate in ("chr" + bare_name, bare_name):
        if candidate in refs:
            return candidate
    return None

# ADDED by Damian, made with ChatGPT. Edited by Teresa with help of ClaudeAI (09-11-26).
def count_and_filter_reads(bam_path, chrom, max_mismatch, min_read_length):
    """
    Count and filter reads per chromosome in a BAM file, tracking discarded reads based on max_mismatch and min_read_length.

    Args:
        bam_path (str): Path to the BAM file.
        chrom (str): Chromosome number.
        max_mismatch (int): Maximum allowed mismatches for reads.
        min_read_length (int): Minimum length of reads.

    Returns:
        tuple: Total reads, reads skipped due to being unmapped, discarded reads due to mismatch, secondary/supplementary reads skipped due to None type, discarded reads due to length, and retained reads.
    """
    # bamfile = pysam.AlignmentFile(bam_path, "rb")

    chrom_label = "chr" + str(chrom).removeprefix("chr")
    total_reads = 0
    skipped_unmapped = 0
    discarded_mismatch = 0
    skipped_Nonetype = 0
    discarded_length = 0
    retained_reads = 0

    # for read in bamfile.fetch():
    #     total_reads += 1

    #     # Skip unmapped reads — they carry no alignment position, no NM tag,
    #     # and query_sequence is None for mate-rescued unmapped entries.
    #     if read.is_unmapped:
    #         continue

    #     # Check mismatch
    #     mismatch_count = read.get_tag('NM') if read.has_tag('NM') else 0
    #     if mismatch_count > max_mismatch:
    #         discarded_mismatch += 1
    #         continue  # Discard read if mismatch count exceeds max_mismatch

    #     # query_sequence is None for secondary/supplementary alignments that
    #     # omit the sequence field; skip them rather than crashing on len().
    #     if read.query_sequence is None:
    #         continue

    #     # Check read length
    #     if len(read.query_sequence) < min_read_length:
    #         discarded_length += 1
    #         continue  # Discard read if it is shorter than the minimum length

    #     retained_reads += 1

    # bamfile.close()

    # Altered version of filtering step above that incorporates resolve_contig() (Teresa -- 09-11-26).
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contig = resolve_contig(bam, chrom)
        # Ensure contig for chromosome in the BAM file is correct.
        if contig is None:
            logger.warning(f"Contig for {chrom_label} not found in {bam_path}.")
            return None
        for read in bam.fetch(contig):
            total_reads += 1

            # Skip unmapped reads — they carry no alignment position, no NM tag, and query_sequence is None for mate-rescued unmapped entries.
            if read.is_unmapped:
                skipped_unmapped += 1
                continue

            # Check mismatch: NM first, then nM (STARsolo/CellRanger), as in BamFilter() in germline.py.
            # Reads without either tag cannot be judged on mismatches and are not discarded here.
            if read.has_tag('NM'):
                mismatch_count = read.get_tag('NM')
            elif read.has_tag('nM'):
                mismatch_count = read.get_tag('nM')
            else:
                mismatch_count = None
            if mismatch_count is not None and mismatch_count >= max_mismatch:
                discarded_mismatch += 1
                continue  # Discard read if mismatch count meets or exceeds max_mismatch

            # query_sequence is None for secondary/supplementary alignments that omit the sequence field; skip them rather than crashing on len().
            if read.query_sequence is None:
                skipped_Nonetype += 1
                continue

            # Check read length
            if len(read.query_sequence) < min_read_length:
                discarded_length += 1
                continue  # Discard read if it is shorter than the minimum length

            retained_reads += 1

    # Log the filtering results # Combined seperate logger.info() lines into one to avoid results from different chromosomes overlapping (Teresa 09-11-26).
    logger.info(
        f"BAM File: {bam_path}\n"
        f"Chromosome: {chrom_label}\n"
        f"Total Reads: {total_reads}\n"
        f"Discarded unmapped sequences: {skipped_unmapped}\n"
        f"Discarded (mismatch >= {max_mismatch}): {discarded_mismatch}\n"
        f"Discarded 'Nonetype' sequences: {skipped_Nonetype}\n"
        f"Discarded (length < {min_read_length}): {discarded_length}\n"
        f"Retained Reads: {retained_reads}"
    )

    return total_reads, skipped_unmapped, discarded_mismatch, skipped_Nonetype, discarded_length, retained_reads

# Function to perform pre-processing of bam files -- 2024-08-08
def preProcess(args):
    """
    Perform pre-processing of BAM files before variant calling.
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """

    logger.info("Performing data preprocess before variant calling with the following arguments.")
    print_parameters_given(args)

    # UMI collapsing needs a UMI tag and a cell barcode tag in the BAM; SmartSeq2 has no UMIs and
    # CEL-Seq2 pipelines keep the UMI in the read name, so stop before doing any work -- 2026-11-09
    if args.umi_collapse and args.platform_library != "10x":
        print(f"ERROR: --umi-collapse is only supported with --platform-library 10x (got '{args.platform_library}'). "
              f"SmartSeq2 has no UMIs and CEL-Seq2 stores the UMI in the read name, not in a BAM tag.")
        sys.exit(1)

    if not os.path.isfile(args.bamFile):
        print(f"ERROR: The list of bam file(s) '{args.bamFile}' cannot be found!")
        sys.exit(1)

    # Create necessary directories -- 2024-09-17
    if args.verbose:
        print(f"\n> Checking the existence of the necessary output directories. If they do not exist, they will be created.")
    if not args.out:
        print(f"ERROR: Output directory not specified!")
        sys.exit(1)
    prepare_output_dirs(args.out, ['preprocess', 'scripts'], verbose=args.verbose)

    sample = []
    # Check the existence of the bam files -- 2024-09-17
    if args.verbose:
        print(f"> Checking the existence of the bam files.")
    with open(args.bamFile) as f_in:
        for line in f_in:
            record = line.strip().split(",")
            sample.append(record[0])
            if args.verbose:
                print(f"> Checking sample {record[0]}")
            logger.debug("Checking sample {}".format(record[0]))
            if len(record) != 2:
                print(f"ERROR: Line with sample name '{record[0]}' does not have exactly 2 comma-delimited columns!")
                sys.exit(1)
            bam_path = Path(record[1])
            if not bam_path.is_file():
                print(f"ERROR: BAM file '{record[1]}' cannot be found!")
                sys.exit(1)
            if not (bam_path.parent / (bam_path.name + ".bai")).is_file() and not (bam_path.with_suffix(bam_path.suffix + ".bai")).is_file():
                # Try both conventions: file.bam.bai and file.bai
                print(f"ERROR: BAM index '{record[1]}.bai' is missing!")
                sys.exit(1)

    # Handle UMI collapse and UMI tag -- 2024-09-18
    if args.umi_collapse:
        umi_collapse = True
        umi_tag = args.umi_collapse if isinstance(args.umi_collapse, str) else "UB"
    else:
        umi_collapse = False
        umi_tag = None

    para_lst = []
    # Process each sample in parallel -- 2024-09-18
    # A list of bam files for each sample is expected
    # sample1, /path/to/sample1.bam
    # sample2, /path/to/sample2.bam
    with open(args.bamFile) as f_in:
        for line in f_in:
            record = line.strip().split(",")
            if args.verbose:
                print(f"> PreProcessing sample {record[0]}")
            logger.debug("PreProcessing sample {}".format(record[0]))
            # Process chromosome 1-22 -- 2024-09-16
            if args.verbose:
                print(f"  - Processing chromosomes 1-22, X...")
            for chrom in range(1, 23):
                if args.debug:
                    print(f"  -- DEBUGGING: chromosome [{chrom}]")
                para_single = dict(
                    chr="chr" + str(chrom), # Add the chromosome number to the chromosome name, make strings of the numbers
                    out=str(Path(args.out)),
                    id=record[0], # record[0] is the sample ID and this assigned here
                    bamFile=record[1], # record[1] is the path to the bam file and this assigned here
                    max_mismatch=args.max_mismatch,
                    samtools=samtools,
                    platform_library=args.platform_library,  # new code -- 2024-09-16
                    verbose=args.verbose,  # new code -- 2024-09-17
                    debug=args.debug,  # new code -- 2024-09-17
                    min_read_length=args.min_read_length,  # new code -- 2024-09-18
                    umi_collapse=umi_collapse,  # Pass umi_collapse flag
                    umi_tag=umi_tag  # Pass umi_tag if provided
                )
                para_lst.append(para_single)

                # # --- Quick read count --- ADDED by Damian, made with ChatGPT.
                # bam_path = record[1]
                # try:
                #     bamfile = pysam.AlignmentFile(bam_path, "rb")
                #     read_count = bamfile.count(
                #         contig="chr" + str(chrom))
                #     bamfile.close()
                #     print(
                #         f"BAM {bam_path} has {read_count} reads on chromosome chr{chrom}")
                #     # Here we call count_and_filter_reads for each BAM file, which already logs stats
                #     count_and_filter_reads(bam_path, args.max_mismatch,
                #                            args.min_read_length)
                # except ValueError as e:
                #     print(
                #         f"Could not count reads for chr{chrom} in BAM {bam_path}: {e}")

            # Process chromosome X -- 2024-09-17
            for chrom in ["X"]:
                if args.debug:
                    print(f"  -- DEBUGGING: chromosome [{chrom}]")
                para_single = dict(
                    chr="chr" + chrom, # Add the chromosome number to the chromosome name
                    out=str(Path(args.out)),
                    id=record[0], # record[0] is the sample ID and this assigned here
                    bamFile=record[1], # record[1] is the path to the bam file and this assigned here
                    max_mismatch=args.max_mismatch,
                    samtools=samtools,
                    platform_library=args.platform_library,  # new code -- 2024-09-16
                    verbose=args.verbose,  # new code -- 2024-09-17
                    debug=args.debug,  # new code -- 2024-09-17
                    min_read_length=args.min_read_length,  # new code -- 2024-09-18
                    umi_collapse=umi_collapse,  # Pass umi_collapse flag
                    umi_tag=umi_tag  # Pass umi_tag if provided
                )
                para_lst.append(para_single)

                # # --- Quick read count --- ADDED by Damian, made with ChatGPT.
                # bam_path = record[1]
                # try:
                #     bamfile = pysam.AlignmentFile(bam_path, "rb")
                #     read_count = bamfile.count(contig="chr" + chrom)
                #     bamfile.close()
                #     print(
                #         f"BAM {bam_path} has {read_count} reads on chromosome chr{chrom}")
                #     # Here we call count_and_filter_reads for each BAM file, which already logs stats
                #     count_and_filter_reads(bam_path, args.max_mismatch,
                #                            args.min_read_length)
                # except ValueError as e:
                #     print(
                #         f"Could not count reads for chr{chrom} in BAM {bam_path}: {e}")

    # Run the BamFilter function in parallel for each chromosome from the germline.py module -- 2024-09-16
    with Pool(processes=args.nthreads) as pool:
        result = pool.map(BamFilter, para_lst)

    # Relocating read count based on new count_and_filter_reads() function (Teresa with help of ClaudeAI -- 09-11-2026)
    for para in para_lst:
        count_and_filter_reads(para["bamFile"], para["chr"], para["max_mismatch"], para["min_read_length"])

    # NEW CODE: Create bam file list for chromosomes 1-22, X
    if args.verbose:
        print(f"> Creating the bam file list for each chromosome, 1-22, X.")
    out_path = Path(args.out)
    for chr in list(range(1, 23)) + ["X"]:
        if args.debug:
            print(f"  -- DEBUGGING: chromosome [{chr}]")
        bamlist_path = out_path / "Bam" / f"chr{chr}.filter.bam.lst"
        with bamlist_path.open("w") as bamlist:
            for s in sample:
                bam_file = out_path / "Bam" / f"{s}_chr{chr}.filter.bam"
                bamlist.write(str(bam_file) + "\n")

    # Disabled -- 2026-11-09: args.study never exists, so this wrote a runGermline_study.sh that sorted a non-existent VCF.
    # # Write a SLURM/script file for preprocessing using the shared helper
    # jobid = args.study if hasattr(args, 'study') else 'study'
    # command = f'echo "[MonopogenLite.py --preProcess] Done preprocessing {jobid}."'
    # write_job_script(jobid, command, out_path, version=VERSION, verbose=args.verbose)

# Main function -- 2024-08-08
def main():
    """
    Main function to parse command-line arguments and execute
    the appropriate subcommand.
    """

    parser = argparse.ArgumentParser(
        description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
MonopogenLite: SNV calling and phasing from single-cell sequencing data.""",
        epilog=f"""\n
Available subcommands: preProcess, germline\n
To see help for each subcommand, use:\n
python MonopogenLite.py preProcess --help\n
python MonopogenLite.py germline --help\n\n
+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} + \n
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter)

    # Add --logfile argument to the main parser BEFORE subparsers
    parser.add_argument('--logfile', required=False,
                        help="Optional log file path. If not provided, logs will be saved to OUTDIR/MonopogenLite.YYYYMMDD.log")

    # Create a subparser object for subcommands (e.g., preProcess, germline)
    subparsers = parser.add_subparsers(title='Subcommands', dest="subcommand", help="Choose one of the available subcommands.")
    # Define a common parser for arguments that both subcommands share (if needed)
    common_parser = argparse.ArgumentParser(add_help=False)

    # Add the subcommand for 'preProcess'
    parser_preProcess = subparsers.add_parser('preProcess', parents=[common_parser],
    help='Preprocess BAM files before variant calling.',
    description='Preprocess of BAM files including removing reads with high alignment mismatches. Default mismatch threshold is 3.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_preProcess.add_argument('-b', '--bamFile', required=True,
                                    help="The comma-separated listf of bam-files for the study sample. The first column should have the sampleID, and the second column the location of the corresponding bam-file. The bam-files should be sorted and indexed. If there are multiple samples, each row with each sample. Required.")
    parser_preProcess.add_argument('-o', '--out', required=True,
                                    help="The output directory. The output will be saved in the output directory. Required.")
    parser_preProcess.add_argument('-a', '--app-path', required=True,
                                    help="The app library paths used in the tool. The app library paths should include (a symlink to) beagle. Also see wiki for installation instructions of relevant tools (samtools, bcftools, bgzip, and java). Required.")
    parser_preProcess.add_argument('-m', '--max-mismatch', required=False, type=int, default=3,
                                    help=(
                                        "Mismatch threshold for reads used in variant calling: reads with this many or more "
                                        "mismatches (NM tag, or nM for STAR-based BAMs) are removed. "
                                        "The default of 3 keeps reads with at most 2 mismatches, as in the original Monopogen."
                                    ))
    parser_preProcess.add_argument('-t', '--nthreads', required=False, type=int, default=1,
                                    help="Number of threads used for SNVs calling. Default is 1.")
    parser_preProcess.add_argument('-l', '--platform-library', required=True, choices=['10x','smartseq2','celseq2'],
                                    help=(
                                        "The platform library used for sequencing. Required. "
                                        "Use '10x' for 10x Genomics RNA-seq (CellRanger) and 10x ATAC-seq (CellRanger ATAC). "
                                        "Use 'smartseq2' for SmartSeq2/STARsolo full-length RNA-seq. "
                                        "Use 'celseq2' for CEL-Seq2 (e.g. bowtie2-based pipelines)."
                                    ))
    parser_preProcess.add_argument('-r', '--min-read-length', required=False, type=int, default=30,
                                    help="The minimum read length for variant calling. Default is 30.")
    parser_preProcess.add_argument('-u', '--umi-collapse', required=False, nargs='?', const=True, default=None,
                                    help=(
                                        "Collapse PCR duplicates by UMI. Not recommended for germline variant calling: "
                                        "duplicate reads add depth, and base and mapping quality filters in bcftools handle PCR errors. "
                                        "Only allowed with --platform-library 10x, for 10x RNA-seq (CellRanger). "
                                        "Omit this flag for: 10x ATAC-seq (no UMIs), SmartSeq2 (no UMIs), and CEL-Seq2 "
                                        "(UMI is encoded in the read name as 'BC-<barcode>_UMI-<umi>', not a BAM tag). "
                                        "Reads are treated as duplicates when they share the cell barcode (CB tag), UMI, "
                                        "start position and strand; the read with the highest mapping quality is kept. "
                                        "Reads without a UMI or cell barcode tag are kept as they are. "
                                        "Optionally give the BAM tag that holds the corrected UMI; the default is 'UB'."
                                    ))
    parser_preProcess.add_argument('-v', '--verbose', action='store_true',
                                    help="Increase output verbosity.")
    parser_preProcess.add_argument('-d', '--debug', action='store_true',
                                    help="For debugging, specifically for installed tools.")
    parser_preProcess.set_defaults(func=preProcess)

    # Set the default function to run when 'preProcess' is called
    parser_preProcess.set_defaults(func=preProcess)

    # Add the subcommand for 'germline'
    parser_germline = subparsers.add_parser('germline', parents=[common_parser],
    help='Germline variant genotype calling and phasing.',
    description='Perform germline variant calling and phasing from single-cell sequencing data.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_germline.add_argument('-r', '--region', required= True,
                                help="The genome regions for variant calling. This file should have either 1 column (chromosome) or 3 columns (chromosome, start, end), where chromosome X is noted as chrX. Required.")
    parser_germline.add_argument('-s', '--step', required=False, default="all", choices=['varScan', 'varProb' , 'varPhasing', 'all'],
                                help=("Run germline variant calling step by step. "
                                      "varScan: call variants with bcftools (and sort and index the unphased VCF); "
                                      "varProb: genotype probabilities with Beagle (needs the varScan output); "
                                      "varPhasing: phase genotypes with Beagle (needs the varProb output); "
                                      "all: all three steps. Default is all."))
    parser_germline.add_argument('--preprocess-dir', required=False, default=None,
                                help=("Directory holding one preProcess output folder per sample (each with a Bam/ folder). "
                                      "Used to build OUT/ListedBamsPerChr/chrN.lst when that list does not exist yet. "
                                      "Default: the --out directory."))
    parser_germline.add_argument('-o', '--out', required= False,
                                help="The output directory. The output will be saved in the output directory. Required.")
    parser_germline.add_argument('-g', '--reference', required= True,
                                help="The human genome reference used for alignment")
    parser_germline.add_argument('-p', '--imputation-panel', required= True,
                                help="The population-level variant panel for variant phasing, such as 1000 Genome phase 3 high-coverage b38 data.")
    parser_germline.add_argument('-m', '--max-softClipped', required=False, type=int, default=1,
                                help="The maximal soft-clipped allowed in one reads for variant calling")
    parser_germline.add_argument('-a', '--app-path', required=True,
                                help="The app library paths used in the tool")
    parser_germline.add_argument('-t', '--nthreads', required=False, type=int, default=1,
                                help="Number of jobs used for SNVs calling")
    parser_germline.add_argument('-n', '--norun', action='store_true',
                                help="Generate the job scripts only. The jobs will not be run.")
    parser_germline.add_argument('-v', '--verbose', action='store_true',
                                help="Increase output verbosity")
    parser_germline.add_argument('-d', '--debug', action='store_true',
                                help="For debugging, specifically for installed tools and some intermediate steps.")
    parser_germline.add_argument('-c', '--chromosome', required=True,
                                 help="Chromosome to call: 1-22 or X, with or without 'chr'. All regions of this chromosome in --region are processed. Required.")

    # Set the default function to run when 'germline' is called
    parser_germline.set_defaults(func=germline)

    # Add --version argument to the main parser
    parser.add_argument('--version', action='version', version=f'{VERSION_NAME} {VERSION} ({VERSION_DATE})')

    # Parse the arguments
    args = parser.parse_args()

    # If no subcommand is provided, print the help and exit
    if args.subcommand is None:
        # if no command is specified, print help and exit
        print("Please specify one subcommand ('preProcess' or 'germline')! Exiting!")
        print("-"*80)
        parser.print_help()
        sys.exit(1)

    # Setup file logging (before subcommand logic)
    today_str = datetime.today().strftime("%Y%m%d")
    resolved_out = str(Path(args.out).resolve() if hasattr(args, 'out') and args.out else Path(".").resolve())
    default_logfile = str(Path(resolved_out) / f"MonopogenLite.{today_str}.log")
    logfile_path = args.logfile or default_logfile
    # Ensure log directory exists
    if logfile_path:
        Path(logfile_path).parent.mkdir(parents=True, exist_ok=True)

    setup_logger(logfile_path)
    if getattr(args, "verbose", False):
        print(f"> Log file set to: {logfile_path}")

    # Log tool versions if verbose
    if getattr(args, "verbose", False):
        log_tool_versions(logging.getLogger())

    # Updated code -- 2024-09-17
    # set paths for tools
    global samtools, bcftools, vcftools, bgzip, java, beagle

    # Execute the shell command to find the location of samtools, bcftools, vcftools, bgzip, and java
    print(f"Checking the existence of the necessary tools.")
    # samtools
    try:
        location_samtools = subprocess.check_output(['which', 'samtools']).strip().decode('utf-8')
        # If you're on Windows, you may need to use where command instead of which.
        # location = subprocess.check_output(['where', 'samtools']).strip().decode('utf-8')
        samtools = os.path.abspath(location_samtools)
        if args.verbose:
            print(f"> samtools location:", samtools)
        if args.debug:
            print(f"> samtools version:", subprocess.check_output([samtools, '--version']).strip().decode('utf-8'))
    except subprocess.CalledProcessError:
        print("ERROR: [samtools] not found.")
    # bcftools
    try:
        location_bcftools = subprocess.check_output(['which', 'bcftools']).strip().decode('utf-8')
        bcftools = os.path.abspath(location_bcftools)
        if args.verbose:
            print(f"> bcftools location:", bcftools)
        if args.debug:
            print(f"> bcftools version:", subprocess.check_output([bcftools, '--version']).strip().decode('utf-8'))
    except subprocess.CalledProcessError:
        print("ERROR: [bcftools] not found.")
    # vcftools
    try:
        location_vcftools = subprocess.check_output(['which', 'vcftools']).strip().decode('utf-8')
        vcftools = os.path.abspath(location_vcftools)
        if args.verbose:
            print(f"> vcftools location:", vcftools)
        if args.debug:
            print(f"> vcftools version:", subprocess.check_output([vcftools, '--version']).strip().decode('utf-8'))
    except subprocess.CalledProcessError:
        print("ERROR: [vcftools] not found.")
    # bgzip
    try:
        location_bgzip = subprocess.check_output(['which', 'bgzip']).strip().decode('utf-8')
        bgzip = os.path.abspath(location_bgzip)
        if args.verbose:
            print(f"> bgzip location:", bgzip)
        if args.debug:
            print(f"> bgzip version:", subprocess.check_output([bgzip, '--version']).strip().decode('utf-8'))
    except subprocess.CalledProcessError:
        print("ERROR: [bgzip] not found.")
    # java
    try:
        location_java = subprocess.check_output(['which', 'java']).strip().decode('utf-8')
        java = os.path.abspath(location_java)
        if args.verbose:
            print(f"> java location:", java)
        if args.debug:
            print(f"> java version:", subprocess.check_output([java, '-version']).strip().decode('utf-8'))
    except subprocess.CalledProcessError:
        print("ERROR: [java] not found.")

    # beagle
    location_beagle = str(Path(args.app_path).resolve() / "beagle.jar")
    beagle = str(Path(location_beagle).resolve())
    if args.verbose:
        print(f"> beagle location:", location_beagle)
    if not Path(beagle).is_file():
        print(f"ERROR: Beagle jar not found at expected path: {beagle}")
        sys.exit(1)

    # Execute the corresponding function for the subcommand
    args.func(args)

    logger.info("Success! See instructions above.")

# Main function -- 2024-08-08
if __name__ == "__main__":
    main()
