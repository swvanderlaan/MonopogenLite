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
VERSION = '1.3.6'
VERSION_DATE = '2025-06-11'
COPYRIGHT = 'Copyright 1979-2025. Jinzhuang Dou | jdou1 [at] mdanderson [dot] org; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
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
def build_sample_commands(samples, chr, out_path, args):
	commands = []
	for sample in samples:
		bam_path = sample["bam"]
		sampleID = sample["sampleID"]
		output_vcf = out_path / "VCF" / f"{sampleID}.chr{chr}.bcf"
		bcftools_cmd = generate_bcftools_command(bam_path, sampleID, chr, out_path)
		if hasattr(args, 'impute') and args.impute:
			if args.genotype:
				beagle_cmd = generate_beagle_cmd_gt(output_vcf, args.app_path, out_path)
			else:
				beagle_cmd = generate_beagle_cmd_gp(output_vcf, args.app_path, out_path)
			command = f"{bcftools_cmd} && {beagle_cmd}"
		else:
			command = bcftools_cmd
		commands.append((sampleID, command))
	return commands

# Function to write a job script for germline variant calling
def write_job_script(jobid, command, out_path, version="1.2.7", verbose=False):
	script_path = out_path / "scripts" / f"runGermline_{jobid}.sh"
	slurm_script_content = f"""#!/bin/bash
echo "[MonopogenLite germline.py v{version}] Start time: $(date)"
echo "[MonopogenLite germline.py v{version}] Running job: {jobid}"

source ~/.bashrc
conda activate monopogen

{command}

echo "[MonopogenLite germline.py v{version}] End time: $(date)"
"""
	write_slurm_script(script_path, slurm_script_content, verbose=verbose)

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
	if verbose:
		print(f"  - SLURM script written to: {path}")

# Function to generate bcftools command for mpileup and filtering
def generate_bcftools_command(bam_filter, jobid, reference, out):
	"""
	Generate the bcftools command for mpileup and filtering.
	Args:
		bam_filter (str): Path to the list of BAM files.
		jobid (str): Job ID or region to process.
		reference (str): Path to the reference genome FASTA file.
		out (str): Output directory where the VCF file will be saved.
	"""

	return (
		f"{bcftools} mpileup -b {bam_filter} --fasta-ref {reference} --regions {jobid} --min-MQ 20 --min-BQ 20 --annotate FORMAT/DP "
		f"| {bcftools} norm -m-both --rm-dup both --check-ref wx --fasta-ref {reference} "
		f"| {bcftools} annotate --set-id '%CHROM:%POS:%REF:%ALT' "
		f"| {bcftools} filter --exclude 'ALT !~ \"^[ATGC]$\"' "
		f"| {bcftools} +fill-tags -Oz -o {out}/germline/{jobid}.gl.vcf.gz"
	)

# Function to generate Beagle command for genotype probabilities
def generate_beagle_cmd_gp(jobid, out, imputation_vcf, chrom, nthreads_downsample):
	"""
	Generate the Beagle command for genotype probabilities.
	Args:
		jobid (str): Job ID or region to process.
		out (str): Output directory where the genotype probabilities will be saved.
		imputation_vcf (str): Path to the imputation panel VCF file.
		chrom (str): Chromosome to process.
		nthreads_downsample (int): Number of threads to use for downsampling.
	"""

	return (
		f"{java} -Xmx20g -jar {beagle} gl={out}/germline/{jobid}.gl.vcf.gz ref={imputation_vcf} "
		f"chrom={chrom} out={out}/germline/{jobid}.gp impute=false modelscale=2 "
		f"nthreads={nthreads_downsample} gprobs=true niterations=0"
	)

# Function to generate Beagle command for genotype phasing
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
		f"{java} -Xmx20g -jar {beagle} gt={out}/germline/{jobid}.gp.vcf.gz ref={imputation_vcf} "
		f"chrom={chrom} out={out}/germline/{jobid}.phased impute=false modelscale=2 "
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
                bam_file_path = Path(out) / "Bam" / f"{chrom}.filter.bam.lst"
                with bam_file_path.open() as bf:
                    for line in bf:
                        bam = line.strip()
                        sample_set.add(bam)

    samples = [{"bam": bam} for bam in sorted(sample_set)]
    return samples

# Function to perform germline variant calling -- 2024-08-08
def germline(args):
	"""
	Perform germline variant calling and phasing from single-cell sequencing data.
	Args:
		args (argparse.Namespace): Parsed command-line arguments.
	"""

	logger.info("Performing germline variant calling with the following arguments.")
	print_parameters_given(args)

	logger.info("Checking existence of essenstial resource files.")
	validate_user_setting_germline(args)

	logger.info("Checking dependencies.")
	check_dependencies(args)

	# Create necessary directories -- 2024-09-17
	if not args.out:
		print(f"ERROR: Output directory not specified!")
		sys.exit(1)
	prepare_output_dirs(args.out, ['germline', 'scripts'], verbose=args.verbose)
	out_path = Path(args.out)

	# Load samples from the region file
	samples = read_sample_list_file(region_file_path=args.region, out=args.out)
	# Convert to expected sample dict format
	samples = [{"sampleID": Path(bam["bam"]).stem.replace(".filter", "").replace(".bam", ""), "bam": bam["bam"]} for bam in samples]
	if args.debug:
		print(f"  -- DEBUGGING: Sample list loaded with {len(samples)} samples.")
		print(f"  -- DEBUGGING: Sample list: {samples}")
	
	# check whether region files were set correctly 
	if args.verbose:
		print(f"  - Checking the region file [{args.region}].")
	if args.debug:
		print(f"  -- DEBUGGING: Creating job-list for pooled job management.")
	joblst = []
	with open(args.region) as f_in:
		for line in f_in:
			# checking the region file and setting the jobid
			if args.debug:
				print(f"  -- DEBUGGING: Line in the region file: {line}; splitting on ',' ...")
			record = line.strip().split(",")
			if args.verbose:
				print(f"  -- chromosome [{record}].")
			if len(record) == 1:
				jobid = record[0]
				if args.debug:
					print(f"  -- DEBUGGING: The region file has only one column. Setting jobid to the chromosome: [{jobid}].")
			if len(record) == 3:
				jobid = record[0] + ":" + record[1] + "-" + record[2]
				if args.debug:
					print(f"  -- DEBUGGING: The region file has three columns. Setting jobid to the region: [{jobid}].")
			
			# setting the bam file list for the given region
			bam_filter = Path(args.out) / "Bam" / f"{record[0]}.filter.bam.lst"

			# checking number of samples within each given region
			N_sample = 0 
			with bam_filter.open() as p:
				for s in p:
					N_sample = N_sample + 1
			if args.debug:
				print(f"  -- DEBUGGING: The number of samples for the given region (should always be the same for each region): {N_sample}.")

			# NEW code -- 2024-09-17
			# check the imputation panel file
			if args.verbose:
				print(f"  - Checking the imputation panel file: [{args.imputation_panel}].")
			imputation_panel_path = Path(args.imputation_panel)
			if record[0] == "chrX":
				imputation_vcf = imputation_panel_path / f"1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.{record[0]}.vcf.gz"
			elif record[0] in {f"chr{n}" for n in range(1, 23)}:
				imputation_vcf = imputation_panel_path / f"1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.{record[0]}.vcf.gz"
			else: 
				print(f"ERROR: The chromosome {record[0]} is not supported!")
				sys.exit(1)
			if args.verbose:
				print(f"  - Starting germline variant calling for region {jobid}.")
			if args.debug:
				print(f"  -- DEBUGGING: number of threads used: {args.nthreads}, attempt 2-fold downsampling these for phasing.")
			nthreads_downsample=int(args.nthreads/2)

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
			cmd3 = generate_beagle_cmd_gp(jobid, str(out_path), imputation_vcf, record[0], nthreads_downsample)
			
			# NEW code -- 2025-06-05
			if args.verbose:
				print(f"    * germline variant phasing of genotype probabilities")
			if args.debug:
				print(f"    -- DEBUGGING: Command to run: {cmd3}")

			# NEW code -- 2025-06-05
			# germline variant phasing of genotypes
			# cmd5 = java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  jobid + ".germline.vcf"  +  " ref=" +  imputation_vcf  +  "  chrom=" + record[0]  + " out="   +  out + "/germline/" + jobid+ ".phased " + "impute=false modelscale=2 nthreads=" + nthreads_downsample + " gprobs=true niterations=0"
			# cmd5 = cmd5 + "\n" + "rm -v " +  out + "/germline/" +  jobid + ".germline.vcf" 
			cmd5 = generate_beagle_cmd_gt(jobid, str(out_path), imputation_vcf, record[0], nthreads_downsample)
			
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
			commands = build_sample_commands(samples, chr, out_path, args)
			for jobid, (_, command) in enumerate(commands):
				write_job_script(jobid, command, out_path, version=VERSION, verbose=args.verbose)

	# pool the shell scripts in the job list when the norun flag is not set
	if not args.norun:
		with Pool(processes=args.nthreads) as pool:
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

# Function to perform pre-processing of bam files -- 2024-08-08
def preProcess(args):
	"""
	Perform pre-processing of BAM files before variant calling.
	Args:
		args (argparse.Namespace): Parsed command-line arguments.
	"""

	logger.info("Performing data preprocess before variant calling with the following arguments.")
	print_parameters_given(args)

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
		umi_tag = args.umi_collapse if isinstance(args.umi_collapse, str) else "UMI"
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
			for chr in range(1, 23):
				if args.debug:
					print(f"  -- DEBUGGING: chromosome [{chr}]")
				para_single = dict(
					chr="chr" + str(chr), # Add the chromosome number to the chromosome name, make strings of the numbers
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

			# Process chromosome X -- 2024-09-17
			for chr in ["X"]:
				if args.debug:
					print(f"  -- DEBUGGING: chromosome [{chr}]")
				para_single = dict(
					chr="chr" + chr, # Add the chromosome number to the chromosome name
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

	# Run the BamFilter function in parallel for each chromosome from the germline.py module -- 2024-09-16
	with Pool(processes=args.nthreads) as pool:
		result = pool.map(BamFilter, para_lst)

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

	# Write a SLURM/script file for preprocessing using the shared helper
	jobid = args.study if hasattr(args, 'study') else 'study'
	command = f'echo "[MonopogenLite.py --preProcess] Done preprocessing {jobid}."'
	write_job_script(jobid, command, out_path, version=VERSION, verbose=args.verbose)
		
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
									help="The maximal alignment mismatch allowed in one reads for variant calling. Default is 3.")
	parser_preProcess.add_argument('-t', '--nthreads', required=False, type=int, default=1,
									help="Number of threads used for SNVs calling. Default is 1.")
	parser_preProcess.add_argument('-l', '--platform-library', required=True, choices=['10x','smartseq2','celseq2'], 
									help="The platform library used for sequencing. This can be 10x, smartseq2, or celseq2. Required.")	
	parser_preProcess.add_argument('-r', '--min-read-length', required=False, type=int, default=30,
									help="The minimum read length for variant calling. Default is 30.")	
	parser_preProcess.add_argument('-u', '--umi-collapse', required=False, nargs='?', const=True, default=None,
									help="Collapse UMIs. By default no UMIs are collapsed. Optionally, specify the UMI tag (e.g., 'UMI', 'RX' or 'MI'). If no tag is provided, default is 'RX'.")	
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
	parser_germline.add_argument('-s', '--step', required= True, default="all", choices=['varScan', 'varProb' , 'varPhasing', 'all'],
								help="Run germline variant calling step by step. varScan: variant calling; varProb: variant phased genotype probabilities; varPhasing: variant phasing; all: all steps. Default is all.")
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
