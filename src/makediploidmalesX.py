#!/usr/bin/env python3

# Importing required modules
import os
import subprocess
import argparse
import gzip
import logging
from datetime import datetime

# Version information
# Change log:
# * v1.1.2 (2024-09-26): Fixed an issue where the check for changes and reverse options was incorrect.
# * v1.1.1 (2024-09-26): Fixed an issue where the changes-file wasn't created.
# * v1.1.0 (2024-09-25): Fix the issue where the haploid genotypes were removed instead of made diploid. Added the option to reverse changes based on a list of variants and samples. Added a logger to log the changes and reversals.
# * v1.0.0 (2024-09-24): Initial version.
VERSION_NAME = 'MakeDiploidMalesX'
VERSION = '1.1.2'
VERSION_DATE = '2024-09-25'
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

# Set up logging
def setup_logger(script_name, verbose):
    """Setup the logger to log to a file with the date and script name, and also log to the console."""
    # Get current date in the format YYYYMMDD
    date_str = datetime.now().strftime('%Y%m%d')
    log_file = f"{date_str}.{script_name}.log"

    # Initialize the logger
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Create file handler to log to a file
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)

    # Create console handler to print to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Formatter for the log messages
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # Add both handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

# Make diploid males in a VCF file
def modify_vcf(input_vcf, output_vcf, changes_file=None, reverse_file=None, logger=None):
    """Modify male genotypes on chromosome X or reverse the changes based on the input lists."""
    
    logger.info(f"Processing input VCF: {input_vcf}")

    changes = []
    
    if reverse_file:
        logger.info(f"Reversing changes based on: {reverse_file}")
        # Load reverse information from the reverse file
        with gzip.open(reverse_file, 'rt') as rev_file:
            reverse_changes = set(tuple(line.strip().split(',')) for line in rev_file)
    else:
        reverse_changes = None

    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):  # Keep header lines
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            variant_id = fields[2]  # Variant ID
            genotypes = fields[9:]  # Genotype columns start at index 9
            variant = fields[0:9]  # Variant information
            
            for i, genotype in enumerate(genotypes):
                if reverse_file:
                    if (variant_id, str(i)) in reverse_changes:
                        original = genotype[0]
                        genotypes[i] = original if original in ['0', '1'] else genotype
                        logger.debug(f"Reversed genotype at {variant_id} for sample {i}")
                else:
                    if genotype == "0" or genotype == "1":
                        new_genotype = f"{genotype}|{genotype}"
                        genotypes[i] = new_genotype
                        changes.append((variant_id, i))
                        logger.debug(f"Modified {genotype} to {new_genotype} at {variant_id} for sample {i}")
            
            outfile.write('\t'.join(variant + genotypes) + '\n')

    # Write changes if not in reverse mode
    if not reverse_file and changes_file:
        logger.info(f"Saving changes to: {changes_file}")
        try:
            with gzip.open(changes_file, 'wt') as ch_file:
                for change in changes:
                    ch_file.write(f"{change[0]},{change[1]}\n")
            logger.info(f"Changes saved successfully to {changes_file}")
        except Exception as e:
            logger.error(f"Error saving changes to {changes_file}: {e}")
            raise

    logger.info(f"Modifications completed. Output written to: {output_vcf}")

# Index the output VCF file using bgzip and tabix
def index_vcf(output_vcf, logger):
    """Index the output VCF file using bgzip and tabix."""
    logger.info(f"Indexing VCF file: {output_vcf}")

    try:
        subprocess.run(["bgzip", "-f", output_vcf], check=True)
        compressed_vcf = output_vcf + ".gz"
        subprocess.run(["tabix", "-f", "-p", "vcf", compressed_vcf], check=True)
        logger.info(f"Indexing completed for {compressed_vcf}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during bgzip/tabix compression or indexing: {e}")
        exit(1)

# Main function
def main():
    # Get the script name (without .py extension)
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Modify or reverse male genotypes on chromosome X in a given VCF-file.""",
		epilog=f"""
For a given input VCF file (`--input-vcf`), this script will change male haploid genotypes 
(0 or 1) to diploid 0|0 or 1|1, and output the modified VCF (`--output`). It can also 
reverse (`--reverse`) the process based on a provided list of samples and variant IDs (`--changes`).
Note that `--changes` and `--reverse` are mutually exclusive, yet `--changes` is required for the
modification to be saved and the output to be reversible using `--reverse`. Thus, the script can be
used to modify a VCF file and save the changes, or reverse the changes based on a list of variants
and samples.

Example to modify:
    python modifyhemigenox.py --input input.vcf --output modified.vcf --changes changes.txt.gz --verbose

Example to reverse:
    python modifyhemigenox.py --input modified.vcf --output original.vcf --reverse changes.txt.gz --verbose

+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} +
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-i", "--input-file", required=True, help="Path to the input VCF file. Required.")
    parser.add_argument("-o", "--output-file", required=True, help="Path to the output modified VCF file. Required.")
    parser.add_argument("-c", "--changes", help="File to save the list of changes (variantID,sample index). Saved as gzipped file (e.g., changes.txt.gz). Optional.")
    parser.add_argument("-r", "--reverse", help="Gzipped file with list of changes to reverse (variantID,sample index). Optional.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output. Optional.")
    parser.add_argument("-V", "--version", action="version", version=f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})")
    
    args = parser.parse_args()
    
    # Check if the changes and reverse options are used together
    if args.changes and args.reverse:
        parser.error("--changes and --reverse cannot be used together.")

    # Check if the input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file does not exist: {args.input_file}")
        exit(1)
    
    # Check if the reverse changes file exists (only when --reverse is provided)
    if args.reverse and not os.path.exists(args.reverse):
        print(f"Error: Reverse file does not exist: {args.reverse}")
        exit(1)
    
    # List the arguments and version through the logger
    logger = setup_logger(script_name, args.verbose)
    logger.info(f"Arguments: {args}")
    logger.info(f"Version: {VERSION_NAME} v{VERSION} ({VERSION_DATE})")

    # Call the modify function
    modify_vcf(args.input_file, args.output_file, changes_file=args.changes, reverse_file=args.reverse, logger=logger)
    
    # Index the output file using bgzip and tabix
    index_vcf(args.output_file, logger=logger)

# Entry point of the script
if __name__ == "__main__":
    main()

