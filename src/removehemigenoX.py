#!/usr/bin/env python3

import argparse
import gzip
import subprocess
import os

# Version information
# Change log:
# v1.0.0 (2024-09-24): Initial version.
VERSION_NAME = 'RemoveHemiGenoX'
VERSION = '1.0.0'
VERSION_DATE = '2024-09-24'
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

def filter_vcf(input_vcf, output_vcf, verbose=False):
    """Filter the VCF file to exclude rows where more than one genotype is not 0|0, 0|1, or 1|1."""
    
    if verbose:
        print(f"Processing input VCF: {input_vcf}")
    
    with gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r') as infile, \
         gzip.open(output_vcf, 'wt') if output_vcf.endswith('.gz') else open(output_vcf, 'w') as outfile:
        
        for line in infile:
            if line.startswith("#"):  # Keep header lines
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            genotypes = fields[9:]  # Genotype columns start at index 9
            valid = True
            
            for genotype in genotypes:
                if not (genotype.startswith("0|0") or genotype.startswith("0|1") or genotype.startswith("1|1") or \
                        genotype.startswith("0/0") or genotype.startswith("0/1") or genotype.startswith("1/1")):
                    valid = False
                    break
            
            if valid:
                outfile.write(line)

        if verbose:
            print(f"Filtering completed. Output written to: {output_vcf}")

def index_vcf(output_vcf, verbose=False):
    """Index the output VCF file using tabix."""
    if verbose:
        print(f"Indexing VCF file: {output_vcf}")

    try:
        subprocess.run(["tabix", "-f", "-p", "vcf", output_vcf], check=True)
        if verbose:
            print(f"Indexing completed for {output_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"Error during tabix indexing: {e}")
        exit(1)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Remove hemizygous genotypes from a given VCF-file.""",
		epilog=f"""
For a given input VCF file (`--input-vcf`), this script will remove all rows where 
more than one genotype is not 0|0, 0|1, or 1|1. The output will be written to a 
new VCF file (`--output-vcf`), which will then be indexed using tabix.

Example:
    python removehemigenoX.py --input input.vcf.gz --output reference.vcf.gz --verbose

+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} +
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-i", "--input-file", required=True, help="Path to the input VCF file (can be gzipped). Should only be applied to chromosome X. Required.")
    parser.add_argument("-o", "--output-file", required=True, help="Path to the output filtered VCF file (gzipped if .gz extension is provided). Required.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output.")
    parser.add_argument("-V", "--version", action="version", version=f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})")
    
    args = parser.parse_args()
    
    # Call the filter function
    filter_vcf(args.input_file, args.output_file, verbose=args.verbose)
    
    # Index the output file using tabix
    index_vcf(args.output_file, verbose=args.verbose)

if __name__ == "__main__":
    main()