#!/usr/bin/env python3
"""
The main interface of MonopgenLite.
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
import subprocess # run shell commands
import multiprocessing as mp
from multiprocessing import Pool

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

# Change log:
# * v1.2.3, 2024-09-19: Added more information the the runGermline scripts.
# * v1.2.2, 2024-09-18: Updated germline.py to properly account for the 'RG' header given different platforms. 
# * v1.2.1, 2024-09-18: Added a separate script to count overlapping variants between input (derived from bam-files) and output VCF files.
# * v1.2.0, 2024-09-18: Added support to provide a minimum read length for variant calling (defaults to 30), added support for the platform library used for sequencing (10x, smartseq2, celseq2), and added support for collapsing UMIs.
# * v1.1.0, 2024-09-18: Added chromosome X support.
# * v1.0.0, 2024-09-19: Initial version. MonoPogenLite is a light-version fork of Monopogen. It only includes the germline-variant-calling pipeline.
# Version and license information
VERSION_NAME = 'MonopogenLite'
VERSION = '1.2.3'
VERSION_DATE = '2024-09-19'
COPYRIGHT = 'Copyright 1979-2024. Jinzhuang Dou | jdou1 [at] mdanderson [dot] org; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
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
LIB_PATH = os.path.abspath(
	os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

if LIB_PATH not in sys.path:
	sys.path.insert(0, LIB_PATH)

PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)

# Function to print errors if any exist -- 2024-08-08
def error_check(all, output, step):
		job_fail = 0
		for id in all:
			if id not in output:
				logger.error("In "+ step + " step " + id + " failed!")
				job_fail = job_fail + 1

		if job_fail > 0:
			logger.error("Failed! See instructions above.")
			sys.exit(1)

# Function to perform germline variant calling -- 2024-08-08
def germline(args):
	logger.info("Performing germline variant calling with the following arguments.")
	print_parameters_given(args)

	logger.info("Checking existence of essenstial resource files.")
	validate_user_setting_germline(args)

	logger.info("Checking dependencies.")
	check_dependencies(args)

	# Create necessary directories -- 2024-09-17
	if args.verbose:
		print(f"Checking the existence of the necessary output directories. If they do not exist, they will be created.")
	if args.debug:
		print(f"  -- DEBUGGING: Output directory: {args.out}.")
	out = args.out
	if args.out:
		os.makedirs(args.out, exist_ok=True)
		if args.verbose:
			print(f"  - Created output directory: {args.out}")
		os.makedirs(os.path.join(args.out, 'germline'), exist_ok=True)
		if args.verbose:
			print(f"  - Created directory to store germline (known) variant data: {os.path.join(args.out, 'germline')}")
		os.makedirs(os.path.join(args.out, 'scripts'), exist_ok=True)
		if args.verbose:
			print(f"  - Created directory to store scripts: {os.path.join(args.out, 'scripts')}")
	else:
		print(f"ERROR: Output directory not specified!")
		sys.exit(1)

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
			if(len(record)==1):
				jobid = record[0]
				if args.debug:
					print(f"  -- DEBUGGING: The region file has only one column. Setting jobid to the chromosome: [{jobid}].")
			if(len(record)==3):
				jobid = record[0] + ":" + record[1] + "-" + record[2]
				if args.debug:
					print(f"  -- DEBUGGING: The region file has three columns. Setting jobid to the region: [{jobid}].")
			
			# setting the bam file list for the given region
			bam_filter = args.out + "/Bam/" +  record[0] +  ".filter.bam.lst"

			# checking number of samples within each given region
			N_sample = 0 
			with open(bam_filter) as p:
				for s in p:
					N_sample = N_sample + 1
			if args.debug:
				print(f"  -- DEBUGGING: The number of samples for the given region (should always be the same for each region): {N_sample}.")

			# NEW code -- 2024-09-17
			# check the imputation panel file
			if args.verbose:
				print(f"  - Checking the imputation panel file: [{args.imputation_panel}].")
			if record[0] == "chrX":
				# imputation_vcf = args.imputation_panel + "1kGP_high_coverage_Illumina." + record[0] + ".filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
				imputation_vcf = args.imputation_panel + "1kGP_high_coverage_Illumina.SNVonly_poly.norm.filtered_af_5e4." + record[0] + ".vcf.gz"
			elif record[0] in [f"chr{n}" for n in range(1, 23)]:
				# imputation_vcf = args.imputation_panel + "1kGP_high_coverage_Illumina." + record[0] + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
				imputation_vcf = args.imputation_panel + "1kGP_high_coverage_Illumina.SNVonly_poly.norm.filtered_af_5e4." + record[0] + ".vcf.gz"
			else: 
				print(f"ERROR: The chromosome {record[0]} is not supported!")
				sys.exit(1)
			if args.verbose:
				print(f"  - Starting germline variant calling for region {jobid}.")
			if args.debug:
				print(f"  -- DEBUGGING: number of threads used: {args.nthreads}, attempt 2-fold downsampling these for phasing.")
			nthreads_downsample=int(args.nthreads/2)

			# NEW COMMANDS with bcftools -- 2024-09-17
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
			cmd1 = bcftools + " mpileup -b " + bam_filter + " --fasta-ref "  + args.reference  + " --regions " +  jobid + " --min-MQ 20 --min-BQ 20 --annotate FORMAT/DP "
			cmd1 = cmd1 + " | " + bcftools  + " norm -m-both --rm-dup both --check-ref wx --fasta-ref " + args.reference
			cmd1 = cmd1 + " | " + bcftools + " annotate --set-id '%CHROM:%POS:%REF:%ALT' " 
			cmd1 = cmd1 + " | " + bcftools +  " filter --exclude 'ALT !~ \"^[ATGC]$\"' "  
			cmd1 = cmd1 + " | " + bcftools +   " +fill-tags -Oz -o " + args.out + "/germline/" +  jobid + ".gl.vcf.gz" 

			# NEW code -- 2024-09-17
			if args.verbose:
				print(f"    * germline variant calling")
			if args.debug:
				print(f"    -- DEBUGGING: Command to run: {cmd1}")

			# NEW code -- 2024-09-17
			# germline variant phasing of genotype probabilities
			cmd3 = java + " -Xmx20g -jar " + beagle +  " gl=" +  out + "/germline/" +  jobid + ".gl.vcf.gz"  +  " ref=" +  imputation_vcf   + "  chrom=" + record[0] + " out="   +  out + "/germline/" + jobid + ".gp " + "impute=false modelscale=2 nthreads=" + str(nthreads_downsample) + " gprobs=true niterations=0"
			
			# NEW code -- 2024-09-17
			if args.verbose:
				print(f"    * germline variant phasing of genotype probabilities")
			if args.debug:
				print(f"    -- DEBUGGING: Command to run: {cmd3}")

			# NEW code -- 2024-09-17
			# germline variant phasing of genotypes
			# cmd5 = java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  jobid + ".germline.vcf"  +  " ref=" +  imputation_vcf  +  "  chrom=" + record[0]  + " out="   +  out + "/germline/" + jobid+ ".phased " + "impute=false modelscale=2 nthreads=" + nthreads_downsample + " gprobs=true niterations=0"
			# cmd5 = cmd5 + "\n" + "rm -v " +  out + "/germline/" +  jobid + ".germline.vcf" 
			
			cmd5 = java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  jobid + ".gp.vcf.gz"  +  " ref=" +  imputation_vcf  +  "  chrom=" + record[0]  + " out="   +  out + "/germline/" + jobid+ ".phased " + "impute=false modelscale=2 nthreads=" + str(nthreads_downsample) + " gprobs=true niterations=0"
			
			# NEW code -- 2024-09-17
			if args.verbose:
				print(f"    * germline variant phasing of genotypes")
			if args.debug:
				print(f"    -- DEBUGGING: Command to run: {cmd5}")

			# NEW code -- 2024-09-17
			# write the commands to a shell script
			if args.verbose:
				print(f"  - Writing the commands to a shell script.")
			f_out = open(out + "/scripts/runGermline_" +  jobid +  ".sh","w")
			if args.step == "varScan" or args.step == "all":
				# NEW code -- 2024-09-17
				if args.verbose:
					print(f"    * variant calling command")
				# writing some extra information to the shell script
				f_out.write({VERSION_NAME} + "v" + {VERSION} + "(" + {VERSION_DATE} + ") \nMonopogenLite: SNV calling and phasing from single-cell sequencing data.\n")
				f_out.write("\nSettings:\n")
				f_out.write("  - Region.................: " + jobid + "\n")
				f_out.write("  - Reference..............: " + args.reference + "\n")
				f_out.write("  - Imputation panel.......: " + args.imputation_panel + "\n")
				f_out.write("  - Maximum soft-clipped...: " + str(args.max_softClipped) + "\n")
				f_out.write("  - App path...............: " + args.app_path + "\n")
				f_out.write("  - Threads................: " + str(args.nthreads) + "\n")
				f_out.write("  - Output.................: " + out + "\n")
				f_out.write("\n")
				f_out.write("\nVariant genotype calling.\n")
				f_out.write(cmd1 + "\n")
			if args.step == "varProb" or args.step == "all":
				# NEW code -- 2024-09-17
				if args.verbose:
					print(f"    * variant phased genotype probabilities command")
				
				# writing some extra information to the shell script
				f_out.write("\nPhasing variant genotype probabilities.\n")
				f_out.write(cmd3 + "\n")
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
			if args.step == "varPhasing" or args.step == "all":
				# NEW code -- 2024-09-17
				if args.verbose:
					print(f"    * variant genotype phasing command")
				# writing some extra information to the shell script
				f_out.write("\nPhasing hard-called phased genotypes (based on genotype probabilities).\n")
				f_out.write(cmd5 + "\n")
			
			# writing some extra information to the shell script
			f_out.write("\n" + {VERSION_NAME} + " v" + {VERSION} + "." + {COPYRIGHT} + "\n" + {COPYRIGHT_TEXT} + "\n")
			# NEW code -- 2024-08-15
			# append jobs to the job list
			if args.verbose:
				print(f"  - Appending the job " + jobid + " to the job list...")
			joblst.append("bash " + out + "/scripts/runGermline_" +  jobid +  ".sh")
	# close the file
	f_out.close()

	if not args.norun == "TRUE":
		with Pool(processes=args.nthreads) as pool:
			print(f"Running the germline variant calling pipeline in a pool.")
			print(f"\nNumber of threads used: {args.nthreads} (downsample for BEAGLE phasing: {nthreads_downsample}\n")
			print(joblst)
			result = pool.map(runCMD, joblst)

# Function to perform pre-processing of bam files -- 2024-08-08
def preProcess(args):
    logger.info("Performing data preprocess before variant calling with the following arguments.")
    print_parameters_given(args)

    assert os.path.isfile(args.bamFile), "The list of bam file(s) {} cannot be found!".format(args.bamFile)

    # Create necessary directories -- 2024-09-17
    if args.verbose:
        print(f"\n> Checking the existence of the necessary output directories. If they do not exist, they will be created.")
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        if args.verbose:
            print(f"  - Created output directory: {args.out}")
        os.makedirs(os.path.join(args.out, 'Bam'), exist_ok=True)
        if args.verbose:
            print(f"  - Created directory to store filtered [.bam]-files: {os.path.join(args.out, 'Bam')}")
    else:
        print("ERROR: Output directory not specified!")
        sys.exit(1)

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
            assert len(record) == 2, "Every line has to have exactly 2 comma-delimited columns! Line with sample name {} does not satisfy this requirement!".format(record[0])
            assert os.path.isfile(record[1]), "Bam file {} cannot be found!".format(record[1])
            assert os.path.isfile(record[1]+".bai"), "Bam file {} has not been indexed!".format(record[1])

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
                    out=args.out,
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
                    out=args.out,
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
    for chr in list(range(1, 23)) + ["X"]:
        if args.debug:
            print(f"  -- DEBUGGING: chromosome [{chr}]")
        bamlist = open(args.out + "/Bam/chr" + str(chr) + ".filter.bam.lst", "w")
        for s in sample:
            bamlist.write(args.out + "/Bam/" + s + "_chr" + str(chr) + ".filter.bam\n")
        bamlist.close()
		
# Main function -- 2024-08-08
def main():
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
	parser_preProcess.add_argument('-o', '--out', required= False,
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

    # Parse the arguments
	args = parser.parse_args()

    # If no subcommand is provided, print the help and exit
	if args.subcommand is None:
		# if no command is specified, print help and exit
		print("Please specify one subcommand ('preProcess' or 'germline')! Exiting!")
		print("-"*80)
		parser.print_help()
		sys.exit(1)

	# Updated code -- 2024-09-17
	# set paths for tools
	global out, samtools, bcftools, vcftools, bgzip, java, beagle 
	out = os.path.abspath(args.out)
	
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
	location_beagle = os.path.abspath(args.app_path) + "/beagle.jar"
	beagle = os.path.abspath(location_beagle)
	if args.verbose:
		print(f"> beagle location:", location_beagle)
	
    # Execute the corresponding function for the subcommand
	args.func(args)

	logger.info("Success! See instructions above.")

# Main function -- 2024-08-08
if __name__ == "__main__":
	main()
