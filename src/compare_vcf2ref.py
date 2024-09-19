#!/usr/bin/env python3

import pysam
import argparse
import os
import logging

# Change log:
# * v1.0.0 2024-09-19: Initial version.
# Version and license information
VERSION_NAME = 'CompareVCFs'
VERSION = '1.0.0'
VERSION_DATE = '2024-09-18'
COPYRIGHT = 'Copyright 1979-2024. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
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

def parse_vcf(vcf_file, debug=False):
    """
    Parse a VCF file and return a set of variants.
    Each variant is represented as a tuple: (chromosome, position, reference allele, alternate allele)
    """
    variants = set()
    if debug:
        logger.info(f"Reading variants from VCF file: {vcf_file}.")
    with pysam.VariantFile(vcf_file) as vcf:
        if debug:
            logger.info(f"Found {len(vcf.header.samples)} samples in the VCF file.")
        for record in vcf:
            if debug:
                logger.debug(f"Reading variant: {record.chrom}:{record.pos} {record.ref}>{record.alts}")
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = ','.join([str(alt_allele) for alt_allele in record.alts])  # Handle multiallelic sites
            variant_tuple = (chrom, pos, ref, alt)
            variants.add(variant_tuple)
    return variants

def compare_vcfs(input_vcf, reference_vcf, verbose=False, debug=False):
    """
    Compare two VCF files and return overlapping and non-overlapping variants.
    """
    # Parse both VCF files into sets of variants
    if verbose:
        logger.info(f"Reading variants from input VCF file: {input_vcf}.")
    input_variants = parse_vcf(input_vcf, debug=debug)
    if verbose:
        logger.info(f"Reading variants from reference VCF file: {reference_vcf}.")
    reference_variants = parse_vcf(reference_vcf, debug=debug)

    # Find overlapping and non-overlapping variants
    if verbose:
        logger.info(f"Getting overlapping variants.")
    overlapping_variants = input_variants.intersection(reference_variants)
    if verbose:
        logger.info(f"Getting non-overlapping variants.")
    non_overlapping_variants = input_variants.difference(reference_variants)

    return overlapping_variants, non_overlapping_variants

def write_variants_to_file(variants, output_file):
    """
    Write a list of variants to an output file.
    """
    with open(output_file, 'w') as f:
        for variant in sorted(variants):
            f.write(f"{variant[0]}\t{variant[1]}\t{variant[2]}\t{variant[3]}\n")

def main():
    parser = argparse.ArgumentParser(
        description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Compare a VCF-file to a common reference VCF-file.""",
		epilog=f"""
For a given input VCF file (`--input-vcf`), this script will compare it to a given 
reference VCF file (`--reference-vcf`) and list the non-overlapping variants. The 
output file will automatically be named as `input_vcf_filename.non_overlapping_variants.txt`.

Optionally, you can enable verbose output (`--verbose`).

Example:
    python compare_vcf2ref.py --input-vcf input.vcf.gz --reference-vcf reference.vcf.gz --verbose

+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} +
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter)    
    parser.add_argument('-i', '--input-vcf', required=True, help="Path to the input VCF file. Required.")
    parser.add_argument('-r', '--reference-vcf', required=True, help="Path to the reference VCF file. Required.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output. Optional.")
    parser.add_argument('-d', '--debug', action='store_true', help="Enable debug output. Optional.")
    parser.add_argument('-V', '--version', action='version', version=f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})")
    
    args = parser.parse_args()
    
    # Setup logging
    logging_level = logging.DEBUG if args.debug else (logging.INFO if args.verbose else logging.WARNING)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging_level)
    logger = logging.getLogger()

    # Derive the output file name from the input VCF file, strip ".vcf.gz" or ".vcf"
    input_vcf_basename = os.path.basename(args.input_vcf)
    output_file = f"{os.path.splitext(input_vcf_basename)[0].replace('.vcf', '').replace('.gz',
