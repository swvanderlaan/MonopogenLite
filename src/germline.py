#!/usr/bin/env python3
"""
The main interface of germline-module.
"""

# Import the required libraries -- 2024-08-15
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
from collections import defaultdict
import multiprocessing as mp
from multiprocessing import Pool

# import data manipulation and analysis libraries
import pandas as pd # data manipulation and analysis
import numpy as np # numerical computing
import gzip # read and write gzip files

# import libraries to handle BAM and VCF files
import pysam # handle BAM files
from pysam import VariantFile # handle VCF files

# Setting the library paths and other global variables -- 2024-08-15
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

# Function to print the parameters
def print_parameters_given(args):
	logger.info("Parameters in effect:")
	for arg in vars(args):
		if arg=="func": continue
		logger.info("--{} = [{}]".format(arg, vars(args)[arg]))

# Function to validate the input sample list
def validate_sample_list_file(args):
	if args.check_hard_clipped:
		out=os.popen("command -v bioawk").read().strip()
		assert out!="", "Program bioawk cannot be found!"

	assert os.path.isfile(args.sample_list), "Sample index file {} cannot be found!".format(args.sample_list)

	try:
		with open(args.sample_list) as f_in:
			for line in f_in:
				record = line.strip().split("\t")
				logger.debug("Checking sample {}".format(record[0]))
				assert len(record)==3, "Every line has to have exactly 3 tab-delimited columns! Line with sample name {} does not satisify this requiremnt!".format(record[0])
				assert os.path.isfile(record[1]), "Bam file {} cannot be found!".format(record[1])
				assert os.path.isfile(record[1]+".bai"), "Bam file {} has not been indexed!".format(record[1])
				assert os.path.isabs(record[1]), "Please use absolute path for bam file {}!".format(record[1])

				if args.check_hard_clipped:
					logger.debug("Checking existence of hard-clipped reads.")
					cmd = "samtools view {} | bioawk -c sam 'BEGIN {{count=0}} ($cigar ~ /H/)&&(!and($flag,256)) {{count++}} END {{print count}}'".format(record[1])
					logger.debug("Command: "+cmd)
					out=os.popen(cmd).read().strip()
					logger.debug("Results: "+out)
					assert out=="0", "Bam file {} contains hard-clipped reads without proper flag (0x100) set! Please use -M or -Y options of BWA MEM!".format(record[1])

				try:
					float(record[2])
					assert 0.0 <= float(record[2]) and float(record[2]) <= 1.0, "Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2])
				except:
					logger.error("Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2]))
					raise ValueError(f"Contamination rate of sample {record[0]} is invalid: {record[2]}")

	except Exception as e:
		logger.error("There is something wrong with the sample index file. Check the logs for more information.")
		logger.exception(e)
		raise

# Function to validate the input region file
def validate_user_setting_germline(args):
	assert os.path.isfile(args.reference), "The genome reference fasta file {} cannot be found!".format(args.reference)
	assert os.path.isdir(args.imputation_panel), "Filtered genotype file of the imputation reference panel {} cannot be found!".format(args.imputation_panel)
	assert os.path.isfile(args.region), "The region file {} cannot be found!".format(args.region)
	# check whether each bam file available	
	for chr in range(1, 23):
		# bamFile = args.out + "/Bam/chr" +  str(chr) +  ".filter.bam.lst"
		bamFile = os.path.join(args.out, "Bam", f"chr{chr}.filter.bam.lst")
		with open(bamFile) as f_in:
			for line in f_in:
				line = line.strip()
				assert os.path.isfile(line), "The bam file {} cannot be found!".format(line)
				assert os.path.isfile(line + ".bai"), "The bam.bai file {} cannot be found!".format(line)
	# check whether region files were set correctly 
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			assert len(record)==3 or len(record)==1, "Every line has to have exactly 3 comma-delimited columns chr1,1,100000 or chr1 (on the whole chromosome)! Line with region {} does not satisify this requiremnt!".format(line)

# Function to check the dependencies
def check_dependencies(args):
	# NEW code -- 2024-08-15
	# these programs are installed via conda/mamba
	programs_to_check = ("vcftools", "bgzip",  "bcftools", "samtools", "java")

	for prog in programs_to_check:
		# NEW code -- 2024-08-15
		location_prog = subprocess.check_output(['which', prog]).strip().decode('utf-8')
		progr_out = os.popen("command -v {}".format(location_prog)).read()
		if args.debug:
			print(f"DEBUGGING: progr_out = {progr_out}")
		assert progr_out != "", "Program {} cannot be found!".format(prog)
	# NEW code -- 2024-08-15
	# these programs are downloaded via the Monopogen github repository
	jars_to_check = ("beagle.jar", ) # keeping the comma to make it a tuple
	for jar in jars_to_check:
		jar_path = os.path.join(args.app_path, jar)
		if args.debug:
			print(f"DEBUGGING: Checking JAR file at {jar_path}")
		assert os.path.isfile(jar_path), "Java jar file {} cannot be found at path {}!".format(jar, jar_path)

# Function to add the chr prefix to the bam file
def addChr(in_bam, samtools, verbose=False):
	# edit the sequence names for your output header
	prefix = 'chr'
	out_bam=in_bam+"tmp.bam"
	input_bam = pysam.AlignmentFile(in_bam,"rb")
	new_head = input_bam.header.to_dict()
	for seq in new_head['SQ']:
		if not seq['SN'].startswith(prefix):
			seq['SN'] = prefix + seq['SN']
	# create output BAM with newly defined header
	with pysam.AlignmentFile(out_bam, "wb", header=new_head) as outf:
		for read in input_bam.fetch():
			prefixed_chrom = prefix + read.reference_name
			a = pysam.AlignedSegment(outf.header)
			a.query_name = read.query_name
			a.query_sequence = read.query_sequence
			a.reference_name = prefixed_chrom
			a.flag = read.flag
			a.reference_start = read.reference_start
			a.mapping_quality = read.mapping_quality
			a.cigar = read.cigar
			a.next_reference_id = read.next_reference_id
			a.next_reference_start = read.next_reference_start
			a.template_length = read.template_length
			a.query_qualities = read.query_qualities
			a.tags = read.tags
			outf.write(a)
	input_bam.close()
	outf.close()
	os.system(samtools + " index " +  out_bam)
	if verbose:
		os.system(" mv -v" + out_bam + " " + in_bam)
		os.system(" mv -v" + out_bam + ".bai  " + in_bam + ".bai")
	else:
		os.system(" mv " + out_bam + " " + in_bam)
		os.system(" mv " + out_bam + ".bai  " + in_bam + ".bai")

# Function to sort and filter the chromosomes -- 2024-09-16
def BamFilter(myargs):
	bamFile = myargs.get("bamFile") # Add input BAM file
	search_chr = myargs.get("chr") # Add chromosome
	samtools = myargs.get("samtools") # Add samtools path
	chr = search_chr # Add chromosome
	id = myargs.get("id") # Add sample ID
	out = myargs.get("out") # Add output directory
	verbose = myargs.get("verbose") # Add verbose option
	debug = myargs.get("debug") # Add debug option

	# Get the mismatch and platform library information
	if verbose:
		logger.info(f"Filtering reads by mismatch: {myargs.get('max_mismatch')}")
	max_mismatch = myargs.get("max_mismatch")

	# Get the platform library information
	if verbose:
		logger.info(f"Platform library: {myargs.get('platform_library')}; platform (PL) is assumed ILLUMINA.")
	platform_library = myargs.get("platform_library")

	# Add option to filter reads by length
	if verbose:
		logger.info(f"Filtering reads by length: {myargs.get('min_read_length')}")
	min_read_length = myargs.get("min_read_length") 
    
	# Get umi_collapse and UMI tag
	# Check also this reference: https://github.com/single-cell-genetics/cellsnp-lite/issues/121
	# Get whether UMI collapsing is enabled
	umi_collapse = myargs.get("umi_collapse", False)
	# Default UMI tag is 'RX'
	umi_tag = myargs.get("umi_tag", "RX")
	if umi_collapse:
		logger.info(f"UMI collapsing is enabled. Using tag {umi_tag} to extract UMIs.")
	
	os.system("mkdir -p " + out +  "/Bam")
	infile = pysam.AlignmentFile(bamFile,"rb")
	contig_names = infile.references
	cnt=0 
	for contig in contig_names:
		if contig.startswith("chr"):
			if debug:
				logger.info("The contig {} contains the prefix 'chr'.".format(contig))
			cnt=cnt+1
	if cnt==0:
		if debug:
			logger.info("Contigs do not contain the prefix 'chr'. Removing 'chr' from search_chr if present.")
		if search_chr.startswith("chr"):
			search_chr = search_chr[3:]

	# Read the header information from the BAM file -- 2024-09-17
	tp = infile.header.to_dict()
	# Check if the read group (RG) information is available in the header
	# If not, add the read group information to the header
	if 'RG' not in tp:
		if verbose: 
			logger.info("Read group information not found in the header.")
		# Check the read group information 
		if debug:
			print(f"Read group information:  {tp}")
		# To avoid the format issue, we update the RG flag based on sample information
		sampleID = os.path.splitext(os.path.basename(myargs["bamFile"]))[0]

		
		# Update the read group (RG) information with the dynamically generated LB_value  -- 2024-09-19
		# - 10x Genomics: Often, the BAM files produced by 10x Genomics pipelines will already have 
		# appropriate read groups. The PU might reflect the flowcell or lane ID, and the LB would 
		# represent the library constructed during the 10x library preparation process.
		# - Smart-seq2 or CEL-Seq2: These platforms involve different library preparation protocols 
		# and sequencing setups. 
		# Set up the sample ID and flowcell/experiment for LB based on platform
		### --- OLD CODE -- possibly remove ---
		# if platform_library == "10x":
		# 	logger.info(f"Platform is {platform_library} -- setting library identifier (LB) to [0.1].")
		# 	tp1 = [{'SM':sampleID,'ID':sampleID, 'LB':0.1, 'PL':"ILLUMINA", 'PU':sampleID}]
		# elif platform_library == "smartseq2" or platform_library == "celseq2":
		# 	logger.info(f"Platform is {platform_library} -- only setting sample ID (SM) and cell/sample identifier (ID).")
		# 	tp1 = [{'SM':sampleID,'ID':sampleID}]
		# tp.update({'RG': tp1})
		# # this debug produces a lot of output
		# if debug:
		# 	print(f"Read group information: {tp1} (original: {tp})")
		# NEW code -- 2025-06-05
		# Initialize the read group (RG) information
		logger.info(f"Initializing read group (RG) information for sample {sampleID}.")
		# Create a list to hold the read group information
		tp1 = [{'SM': sampleID, 'ID': sampleID}]  # Base RG fields

		if platform_library == "10x":
			logger.info(f"Platform is {platform_library} -- adding LB, PU, and PL fields.")
			tp1[0].update({
				'LB': sampleID,  # use sampleID or similar string-based identifier
				'PU': sampleID,
				'PL': "ILLUMINA"
			})
		elif platform_library in ["smartseq2", "celseq2"]:
			logger.info(f"Platform is {platform_library} -- using minimal RG fields.")
		else:
			logger.warning(f"Unknown platform_library: {platform_library} -- using default RG fields.")

		tp.update({'RG': tp1})
		# this debug produces a lot of output
		if debug:
			print(f"Read group information: {tp1} (original: {tp})")
  
	outfile =  pysam.AlignmentFile( out + "/Bam/" +id + "_" + chr + ".filter.bam", "wb", header=tp)

	# NEW code -- 2024-09-18
	# Dictionary to group reads by UMI and position
	umi_dict = defaultdict(list)

	# Process reads and group them by UMI and position
	for s in infile.fetch(search_chr):
		val = None
		if s.has_tag("NM"):
			if debug:
				logger.info(f"Read {s.query_name} has NM tag.")
			val = s.get_tag("NM")
		elif s.has_tag("nM"):
			if debug:
				logger.info(f"Read {s.query_name} has nM tag.")
			val = s.get_tag("nM")

        # Filter by mismatch and read length
		# if val < max_mismatch and s.query_length >= min_read_length:
		if val is not None and val < max_mismatch and s.query_length >= min_read_length:
			if debug: 
				logger.info(f"Read {s.query_name} has {val} mismatches and length {s.query_length}.")
			if umi_collapse and s.has_tag(umi_tag):  # Use the specified UMI tag
				if debug:
					logger.info(f"Read {s.query_name} has UMI tag {umi_tag}.")
				umi = s.get_tag(umi_tag)  # Extract UMI from the specified tag
				pos = (s.reference_name, s.reference_start)
				umi_dict[(umi, pos)].append(s)  # Group reads by UMI and position
			else:
				if debug:
					logger.info(f"Read {s.query_name} does not have UMI tag; no UMI collapsing is applied.")
				outfile.write(s)

	# If UMI collapsing is enabled, process grouped reads
	if umi_collapse:
		if debug: 
			logger.info(f"UMI collapsing is enabled. Processing reads grouped by UMI and position.")
		for (umi, pos), reads in umi_dict.items():
			# Collapsing strategy: Keep the read with the highest mapping quality
			best_read = max(reads, key=lambda x: x.mapping_quality)
			outfile.write(best_read)  # Write the collapsed read

	# close out the files
	infile.close()
	outfile.close()

	# Index the BAM file -- 2024-09-16
	logger.info(f"Indexing the BAM file [{out}/Bam/{id}_{chr}.filter.bam].")
	result = subprocess.run([samtools, "index", out + "/Bam/" + id + "_" + chr + ".filter.bam"], capture_output=True, text=True)
	if result.returncode != 0:
		logger.error(f"ERROR: Error running samtools index: {result.stderr}.")

	# Adding the chr prefix to the BAM file when the contig does not include this-- 2024-09-17
	if cnt == 0:
		# Add the chr prefix to the BAM file -- 2024-09-17
		addChr(out + "/Bam/" +  id+ "_" + chr+ ".filter.bam", samtools, verbose=myargs.get("verbose", False))
	bamfile = out + "/Bam/" +  id+ "_" + chr + ".filter.bam"
	return(bamfile)

# Function to get the tag and otherwise return an error message
def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return f"ERROR: tag {tag_name} not found -- NotFound"

# Function to extract the command errors if any
def runCMD(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.returncode == 0:
		logger.info(f"Command '{cmd}' successfully run.")
		return cmd  # Return the command that was successfully run
	else:
		logger.error(f"ERROR: Command '{cmd}' failed with error: {result.stderr}")
		return None  # Handle failure case if needed
