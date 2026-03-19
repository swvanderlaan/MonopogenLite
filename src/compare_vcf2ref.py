#!/usr/bin/env python3

import pysam
import argparse
import os
import logging
import subprocess
import gzip

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

# Initialize logger globally
logger = logging.getLogger()

def log_version_info():
    """
    Logs the version information at the start of the script.
    """
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})")
    logger.info(f"{COPYRIGHT}")
    logger.info(f"{COPYRIGHT_TEXT.splitlines()[1]}")  # Shortened copyright text for the log

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

def get_vcf_version(input_vcf):
    """
    Extract the VCF version from the input VCF file.
    Default to VCFv4.2 if not found.
    """
    open_func = gzip.open if input_vcf.endswith('.gz') else open  # Use gzip.open for .gz files
    with open_func(input_vcf, 'rt') as f:  # Open as text mode
        for line in f:
            if line.startswith('##fileformat='):
                return line.strip().split('=')[1]
    return 'VCFv4.2'  # Default VCF version

def check_and_index_vcf(vcf_file, verbose=False):
    """
    Check if the VCF file is indexed (i.e., has a .tbi index file for .vcf.gz or .csi for .vcf).
    If not, index the VCF using bcftools index.
    """
    # Determine the expected index file based on whether it's a compressed VCF
    if vcf_file.endswith('.vcf.gz'):
        index_file = vcf_file + '.tbi'
    else:
        index_file = vcf_file + '.csi'

    if not os.path.exists(index_file):
        if verbose:
            logger.info(f"Index for {vcf_file} not found. Indexing with bcftools...")
        # Run bcftools index to create the correct index (tbi or csi)
        result = subprocess.run(['bcftools', 'index', vcf_file], capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"ERROR: Failed to index {vcf_file}: {result.stderr}")
            raise RuntimeError(f"Failed to index {vcf_file}")
        if verbose:
            logger.info(f"Index created for {vcf_file}.")
    else:
        if verbose:
            logger.info(f"Index found for {vcf_file}. No indexing needed.")

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

def write_vcf(variants, output_file, input_vcf):
    """
    Write the non-overlapping variants to a new VCF file with the proper VCF headers.
    """
    logger.info(f"Writing non-overlapping variants to VCF file: {output_file}")
    
    # Get the VCF version from the input file
    vcf_version = get_vcf_version(input_vcf)
    
    with pysam.VariantFile(input_vcf) as template_vcf:
        # Copy the header from the input VCF
        header = template_vcf.header

        # Add VCF version to the header
        header.add_meta('fileformat', vcf_version)
        
        # Create a new VariantFile with the copied header
        with pysam.VariantFile(output_file, 'w', header=header) as vcf_out:
            for variant in sorted(variants):
                chrom, pos, ref, alt = variant
                
                # Construct a new VCF record
                new_record = vcf_out.new_record(
                    contig=chrom,
                    start=pos-1,  # VCF uses 0-based indexing
                    stop=pos,     # stop is the same as start + len(ref)
                    alleles=(ref, alt),
                    id=".",
                    qual=None,
                    filter=None,
                    info={}
                )
                
                # Write the new record to the VCF file
                vcf_out.write(new_record)
    
    # Compress the VCF file with bgzip
    subprocess.run(['bgzip', output_file], check=True)

    # Index the compressed VCF file with tabix
    subprocess.run(['tabix', '-p', 'vcf', output_file + '.gz'], check=True)

def main():
    parser = argparse.ArgumentParser(
        description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Compare a VCF-file to a common reference VCF-file.""",
		epilog=f"""
For a given input VCF file (`--input-vcf`), this script will compare it to a given 
reference VCF file (`--reference-vcf`) and list the non-overlapping variants. The 
output file will automatically be named as `input_vcf_filename.non_overlapping_variants.vcf.gz` and 
saved in the same directory as the input VCF file.

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

    # Log version information
    log_version_info()

    # Ensure both VCF files are indexed
    check_and_index_vcf(args.input_vcf, verbose=args.verbose)
    check_and_index_vcf(args.reference_vcf, verbose=args.verbose)

    # Derive the output file name from the input VCF file, strip ".vcf.gz" or ".vcf"
    input_vcf_basename = os.path.basename(args.input_vcf)
    input_vcf_dir = os.path.dirname(args.input_vcf)
    output_file = os.path.join(input_vcf_dir, f"{os.path.splitext(input_vcf_basename)[0].replace('.vcf', '')}.non_overlapping_variants.vcf")
    
    # Compare VCF files
    overlapping_variants, non_overlapping_variants = compare_vcfs(
        args.input_vcf, args.reference_vcf, verbose=args.verbose, debug=args.debug
    )
    
    # Write non-overlapping variants to the output VCF file
    write_vcf(non_overlapping_variants, output_file, args.input_vcf)

    # Log final counts
    total_input_variants = len(overlapping_variants) + len(non_overlapping_variants)
    total_reference_variants = len(parse_vcf(args.reference_vcf))  # To calculate totals in the reference
    total_overlapping_variants = len(overlapping_variants)
    total_non_overlapping_variants = len(non_overlapping_variants)

    logger.info(f"Total variants in input VCF: {total_input_variants}")
    logger.info(f"Total variants in reference VCF: {total_reference_variants}")
    logger.info(f"Overlapping variants: {total_overlapping_variants}")
    logger.info(f"Non-overlapping variants: {total_non_overlapping_variants}")
    logger.info(f"Non-overlapping variants written to: {output_file}.gz")

if __name__ == "__main__":
    main()
