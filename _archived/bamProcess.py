#!/usr/bin/env python3
"""
The main interface of bamProcess.py
"""

# Import the required libraries -- 2024-09-16
import argparse
import sys
import os
import logging
import shutil
import glob
import re
import pysam
import time
import subprocess
import pandas as pd
from pysam import VariantFile

LIB_PATH = os.path.abspath(
	os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

if LIB_PATH not in sys.path:
	sys.path.insert(0, LIB_PATH)

PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

#import pipelines
#from pipelines import get_cluster_cfgfile
#from pipelines import PipelineHandler

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)

# Function to add chromosome prefix to the BAM file, and index the BAM file; used in germline.py -- 2024-09-16
def addChr(args):
	# edit the sequence names for your output header
	in_bam = args.bamFile
	prefix = 'chr'
	out_bam=in_bam+"tmp.bam"
	print(in_bam)
	input_bam = pysam.AlignmentFile(in_bam,"rb")
	new_head = input_bam.header.to_dict()
	for seq in new_head['SQ']:
		seq['SN'] = prefix  + seq['SN']
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
	os.system("samtools  index " +  out_bam)
	os.system(" mv -v" + out_bam + " " + in_bam)
	os.system(" mv -v" + out_bam + ".bai  " + in_bam + ".bai")

# Function to sort and filter the BAM file; used in germline.py -- 2024-09-16
# ORIGINAL CODE
# def sort_chr(chr_lst):
# 	# sort chr IDs from 1...22, X, Y, and MT -- 2024-09-16
# 	chr_lst_sort = []

# 	# Sort chromosomes 1 to 22
# 	for i in range(1, 23):
# 		i = str(i)
# 		if  i in chr_lst:
# 			chr_lst_sort.append(i)
# 		i_chr = "chr"+i 
# 		if  i_chr in chr_lst:
# 			chr_lst_sort.append(i_chr)
# 	chr_lst = chr_lst_sort 
# 	return chr_lst

def sort_chr(chr_lst):
    # sort chr IDs from 1...22, X
    chr_lst_sort = []
    
    # Sort chromosomes 1 to 22
    for i in range(1, 23):
        i = str(i)
        if i in chr_lst:
            chr_lst_sort.append(i)
        i_chr = "chr" + i
        if i_chr in chr_lst:
            chr_lst_sort.append(i_chr)
    
    # Add chromosomes X
    if "X" in chr_lst:
        chr_lst_sort.append("X")
    if "chrX" in chr_lst:
        chr_lst_sort.append("chrX")
    
	# Return the sorted chromosome list
    return chr_lst_sort

