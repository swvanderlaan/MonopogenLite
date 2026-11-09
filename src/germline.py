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

# BAM tag holding the corrected cell barcode (10x Genomics CellRanger); used to keep UMI collapsing within cells.
CELL_BARCODE_TAG = "CB"

# Shared logger -- 2026-11-09
# MonopogenLite.py uses the same logger name and attaches the console and --logfile handlers to it, so
# messages from this module (BamFilter, runCMD) also reach the log file. No handlers are added here.
LOGGER_NAME = "MonopogenLite"
logger = logging.getLogger(LOGGER_NAME)

# Stop with a clear error message instead of an assert (asserts give a traceback and are skipped with python -O) -- 2026-11-09
def fail(message):
    logger.error(message)
    sys.exit(1)

# Function to print the parameters
def print_parameters_given(args):
    logger.info("Parameters in effect:")
    for arg in vars(args):
        if arg=="func": continue
        logger.info("--{} = [{}]".format(arg, vars(args)[arg]))
    if getattr(args, "debug", False):
        print(f"  -- DEBUGGING: Germline after print_parameters_given()")

# Disabled -- 2026-11-09: unused; it validates a tab-delimited sample list (--sample_list, --check_hard_clipped)
# that MonopogenLite has no options for. Kept as comments while testing.
# # Function to validate the input sample list
# def validate_sample_list_file(args):
#     if args.check_hard_clipped:
#         out=os.popen("command -v bioawk").read().strip()
#         assert out!="", "Program bioawk cannot be found!"
#
#     assert os.path.isfile(args.sample_list), "Sample index file {} cannot be found!".format(args.sample_list)
#
#     try:
#         with open(args.sample_list) as f_in:
#             for line in f_in:
#                 record = line.strip().split("\t")
#                 logger.debug("Checking sample {}".format(record[0]))
#                 assert len(record)==3, "Every line has to have exactly 3 tab-delimited columns! Line with sample name {} does not satisify this requiremnt!".format(record[0])
#                 assert os.path.isfile(record[1]), "Bam file {} cannot be found!".format(record[1])
#                 assert os.path.isfile(record[1]+".bai"), "Bam file {} has not been indexed!".format(record[1])
#                 assert os.path.isabs(record[1]), "Please use absolute path for bam file {}!".format(record[1])
#
#                 if args.check_hard_clipped:
#                     logger.debug("Checking existence of hard-clipped reads.")
#                     cmd = "samtools view {} | bioawk -c sam 'BEGIN {{count=0}} ($cigar ~ /H/)&&(!and($flag,256)) {{count++}} END {{print count}}'".format(record[1])
#                     logger.debug("Command: "+cmd)
#                     out=os.popen(cmd).read().strip()
#                     logger.debug("Results: "+out)
#                     assert out=="0", "Bam file {} contains hard-clipped reads without proper flag (0x100) set! Please use -M or -Y options of BWA MEM!".format(record[1])
#
#                 try:
#                     float(record[2])
#                     assert 0.0 <= float(record[2]) and float(record[2]) <= 1.0, "Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2])
#                 except:
#                     logger.error("Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2]))
#                     raise ValueError(f"Contamination rate of sample {record[0]} is invalid: {record[2]}")
#
#     except Exception as e:
#         logger.error("There is something wrong with the sample index file. Check the logs for more information.")
#         logger.exception(e)
#         raise

# Function to validate the input region file
def validate_user_setting_germline(args):
    if not os.path.isfile(args.reference):
        fail(f"The genome reference fasta file [{args.reference}] cannot be found!")
    if not os.path.isdir(args.imputation_panel):
        fail(f"The imputation panel directory [{args.imputation_panel}] cannot be found!")
    if not os.path.isfile(args.region):
        fail(f"The region file [{args.region}] cannot be found!")
    # check whether each bam file available (outdated worked for per sample).
    # Implementing a new check is redundant, because of the presence of validateChromosomalBams.sh
    # for chr in range(1, 23):
    #     # bamFile = args.out + "/Bam/chr" +  str(chr) +  ".filter.bam.lst"
    #     bamFile = os.path.join(args.out, "Bam", f"chr{chr}.filter.bam.lst")
    #     with open(bamFile) as f_in:
    #         for line in f_in:
    #             line = line.strip()
    #             assert os.path.isfile(line), "The bam file {} cannot be found!".format(line)
    #             assert os.path.isfile(line + ".bai"), "The bam.bai file {} cannot be found!".format(line)
    # check whether region files were set correctly
    with open(args.region) as f_in:
        for line in f_in:
            # Blank lines and comment lines are skipped, as in read_regions_for_chromosome() -- 2026-11-09
            if not line.strip() or line.strip().startswith("#"):
                continue
            record = line.strip().split(",")
            if len(record) not in (1, 3):
                fail(f"Every line of the region file has to have 1 column (chr1, the whole chromosome) or 3 comma-delimited "
                     f"columns (chr1,1,100000); this line does not: [{line.strip()}]")

    if getattr(args, "debug", False):
        print(f"  -- DEBUGGING: Germline after validate_user_setting_germline()")

# Function to check the dependencies
def check_dependencies(args):
    # NEW code -- 2024-08-15
    # these programs are installed via conda/mamba
    programs_to_check = ("vcftools", "bgzip",  "bcftools", "samtools", "java")

    for prog in programs_to_check:
        # Look the program up on PATH; a missing program gives a clear error instead of a traceback -- 2026-11-09
        location_prog = shutil.which(prog)
        if args.debug:
            print(f"DEBUGGING: {prog} = {location_prog}")
        if location_prog is None:
            fail(f"Program [{prog}] cannot be found! Is the 'monopogen' environment activated?")
    # NEW code -- 2024-08-15
    # these programs are downloaded via the Monopogen github repository
    jars_to_check = ("beagle.jar", ) # keeping the comma to make it a tuple
    for jar in jars_to_check:
        jar_path = os.path.join(args.app_path, jar)
        if args.debug:
            print(f"DEBUGGING: Checking JAR file at {jar_path}")
        if not os.path.isfile(jar_path):
            fail(f"Java jar file [{jar}] cannot be found at path [{jar_path}]!")
    if args.debug:
        print(f"  -- DEBUGGING: Germline after check_dependencies()")

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
    # Index and move with explicit error checks; the old os.system("mv -v" + path) call lacked a space,
    # so with --verbose the renamed BAM silently never replaced the original -- 2026-11-09
    result = subprocess.run([samtools, "index", out_bam], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools index failed for {out_bam}: {result.stderr.strip()}")
    os.replace(out_bam, in_bam)
    os.replace(out_bam + ".bai", in_bam + ".bai")
    if verbose:
        logger.info(f"Added the 'chr' prefix to the contig names in [{in_bam}].")

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
    # Get whether UMI collapsing is enabled; MonopogenLite.py passes the UMI tag ('UB' by default) when it is.
    umi_collapse = myargs.get("umi_collapse", False)
    umi_tag = myargs.get("umi_tag")
    if umi_collapse:
        logger.info(f"UMI collapsing is enabled. Using tag {umi_tag} for UMIs and tag {CELL_BARCODE_TAG} for cell barcodes.")

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

    # Give the filtered BAM exactly one read group, named after the sample ID from the --bamFile list -- 2026-11-09
    # bcftools mpileup treats every distinct SM in the read groups as a separate sample. Upstream read
    # groups cannot be trusted for this: when absent, the old code named the sample after the BAM file
    # (CellRanger calls every BAM possorted_genome_bam.bam, so all samples merged into one), and STARsolo
    # SmartSeq2 BAMs can carry one read group per cell (so every cell became a sample). MonopogenLite
    # genotypes individuals, so one filtered BAM is always one sample.
    if 'RG' in tp:
        logger.info(f"Replacing {len(tp['RG'])} upstream read group(s) with a single read group for sample {id}.")
        if debug:
            print(f"Upstream read group information: {tp['RG']}")
    else:
        logger.info(f"No read group information found; adding a read group for sample {id}.")
    read_group = {'ID': id, 'SM': id}
    if platform_library == "10x":
        read_group.update({'LB': id, 'PU': id, 'PL': "ILLUMINA"})
    elif platform_library not in ["smartseq2", "celseq2"]:
        logger.warning(f"Unknown platform_library: {platform_library} -- using minimal read group fields.")
    tp['RG'] = [read_group]

    outfile =  pysam.AlignmentFile( out + "/Bam/" +id + "_" + chr + ".filter.bam", "wb", header=tp)

    # UMI collapsing -- 2026-11-09
    # A UMI identifies one molecule within one cell, so duplicates are reads with the same cell barcode,
    # UMI, start position and strand. Reads arrive sorted by position, so once the start position moves
    # on, every group at the previous position is complete and is written straight away. This keeps only
    # one position's reads in memory and keeps the output sorted (and therefore indexable).
    umi_groups = defaultdict(list)
    umi_group_pos = None
    umi_stats = {"reads_with_umi": 0, "molecules_written": 0}

    def flush_umi_groups():
        for reads in umi_groups.values():
            # Collapsing strategy: keep the read with the highest mapping quality
            outfile.write(max(reads, key=lambda x: x.mapping_quality))
            umi_stats["molecules_written"] += 1
        umi_groups.clear()

    # Process reads
    for s in infile.fetch(search_chr):
        # Skip unmapped reads — they carry no alignment position and no NM tag
        if s.is_unmapped:
            if debug:
                logger.info(f"Read {s.query_name} is unmapped; skipping.")
            continue

        val = None
        if s.has_tag("NM"):
            if debug:
                logger.info(f"Read {s.query_name} has NM tag.")
            val = s.get_tag("NM")
        elif s.has_tag("nM"):
            if debug:
                logger.info(f"Read {s.query_name} has nM tag.")
            val = s.get_tag("nM")

        # Filter by mismatch and read length.
        # When neither NM nor nM is present (val is None — e.g. STARsolo SmartSeq2 reads,
        # or aligners that omit the tag) the read cannot be filtered by mismatch, so it passes.
        if (val is None or val < max_mismatch) and s.query_length >= min_read_length:
            if debug:
                logger.info(f"Read {s.query_name} has {val} mismatches and length {s.query_length}.")
            # Point every kept read at the single read group of this sample (see header above).
            s.set_tag("RG", id, value_type="Z")
            # Write out finished UMI groups before any read at a later position, so the output stays sorted.
            if umi_groups and s.reference_start != umi_group_pos:
                flush_umi_groups()
            if umi_collapse and umi_tag and s.has_tag(umi_tag) and s.has_tag(CELL_BARCODE_TAG):
                if debug:
                    logger.info(f"Read {s.query_name} has UMI tag {umi_tag} and cell barcode tag {CELL_BARCODE_TAG}.")
                umi_group_pos = s.reference_start
                umi_groups[(s.get_tag(CELL_BARCODE_TAG), s.get_tag(umi_tag), s.is_reverse)].append(s)
                umi_stats["reads_with_umi"] += 1
            else:
                if debug and umi_collapse:
                    logger.info(f"Read {s.query_name} lacks the UMI or cell barcode tag; no UMI collapsing is applied.")
                outfile.write(s)

    # Write the UMI groups of the last position
    flush_umi_groups()
    if umi_collapse:
        if umi_stats["reads_with_umi"] == 0:
            logger.warning(f"UMI collapsing was requested, but no reads on {chr} carried both the {umi_tag} and {CELL_BARCODE_TAG} tags; nothing was collapsed.")
        else:
            logger.info(f"UMI collapsing on {chr}: {umi_stats['reads_with_umi']} reads with a UMI were collapsed into "
                        f"{umi_stats['molecules_written']} molecules ({umi_stats['reads_with_umi'] - umi_stats['molecules_written']} duplicates removed).")

    # close out the files
    infile.close()
    outfile.close()

    # Index the BAM file -- 2024-09-16
    # An unindexed BAM cannot be used by bcftools mpileup --regions, so stop here instead of continuing -- 2026-11-09
    logger.info(f"Indexing the BAM file [{out}/Bam/{id}_{chr}.filter.bam].")
    result = subprocess.run([samtools, "index", out + "/Bam/" + id + "_" + chr + ".filter.bam"], capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"ERROR: Error running samtools index: {result.stderr}.")
        raise RuntimeError(f"samtools index failed for {out}/Bam/{id}_{chr}.filter.bam: {result.stderr.strip()}")

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
