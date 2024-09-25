#!/usr/bin/env python3

# Import required modules
import argparse
import os
import subprocess
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import logging
from datetime import datetime

# Change log:
# * v1.0.0 2024-09-25: Initial version. 
# Version and license information 
VERSION_NAME = "Variant Reference Checker"
VERSION = "1.0.0"
VERSION_DATE = "2024-09-25"
COPYRIGHT='Copyright 1979-2024. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science'
COPYRIGHT_TEXT='''
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

# Run bcftools stats
def run_bcftools_stats(input_vcf, output_stats, verbose=False):
    """Run bcftools stats to generate statistics from the VCF file."""
    command = ['bcftools', 'stats', input_vcf, '-o', output_stats]
    
    logger.info(f"Running bcftools stats on {input_vcf}...")
    subprocess.run(command, check=True)

# Run bcftools plot-vcfstats
def run_bcftools_plot(stats_file, output_prefix, verbose=False):
    """Run bcftools plot-vcfstats to generate plots from the stats file."""
    command = ['bcftools', 'plot-vcfstats', stats_file, '--prefix', output_prefix]
    
    logger.info(f"Running bcftools plot-vcfstats with prefix {output_prefix}...")
    subprocess.run(command, check=True)

# Generate custom plots
def generate_plots(input_vcf, output_dir, verbose=False):
    """Generate custom plots using cmcrameri color maps."""
    logger.info(f"Generating custom plots with cmcrameri color maps in {output_dir}...")

    # Example: Generating a color map plot with cmcrameri
    x = [i for i in range(10)]
    y = [i**2 for i in range(10)]

    plt.scatter(x, y, c=y, cmap=cmc.batlow)
    plt.colorbar(label='Color scale: batlow')
    
    plt.title('Sample Plot')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    # Create the base name by removing the .vcf.gz suffix from the input VCF file
    base_name = os.path.basename(input_vcf).replace(".vcf.gz", "")
    
    # Save as PNG and PDF with appropriate file names
    plt.savefig(os.path.join(output_dir, f"{base_name}_sample_plot.png"))
    plt.savefig(os.path.join(output_dir, f"{base_name}_sample_plot.pdf"))

    logger.info(f"Plots saved as PNG and PDF in {output_dir}.")

# Main function
def main():
    parser = argparse.ArgumentParser(description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Run bcftools stats and generate plots using cmcrameri colors.""",
		epilog=f"""
For a given input VCF file (`--input-vcf`), this script will run bcftools stats and
generate plots using cmcrameri color maps. The output files will be saved in the same
directory as the input VCF file.

Example:
    python variant_ref_checker.py --input input.vcf.gz --verbose

+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} +
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', required=True, help="Path to the input VCF file (vcf.gz format). Required.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output.")
    parser.add_argument('-V', '--version', action='version', version=f'{VERSION_NAME} {VERSION} ({VERSION_DATE})')
    
    args = parser.parse_args()
    
    input_vcf = args.input
    output_dir = os.path.dirname(input_vcf)
    
    # Generate the stats file and plot prefixes
    base_name = os.path.basename(input_vcf).replace(".vcf.gz", "")
    output_stats = os.path.join(output_dir, f'{base_name}_vcf_stats.txt')
    output_prefix = os.path.join(output_dir, base_name)

    # Ensure the input file is a .vcf.gz file
    if not input_vcf.endswith(".vcf.gz"):
        raise ValueError("Input file must be a .vcf.gz file.")
    
    # Set up logger
    logger = setup_logger("variant_ref_checker", args.verbose)
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})")
    logger.info(f"Settings:")
    logger.info(f"Input file..........: {input_vcf}")
    logger.info(f"Output directory....: {output_dir}")
    logger.info(f"Output stats file...: {output_stats}")
    logger.info(f"Verbose.............: {args.verbose}")
    logger.info(f"")
    
    # Run bcftools stats
    run_bcftools_stats(input_vcf, output_stats, args.verbose)
    
    # Run bcftools plot-vcfstats
    run_bcftools_plot(output_stats, output_prefix, args.verbose)
    
    # Generate additional custom plots
    # generate_plots(input_vcf, output_dir, args.verbose)

# Run the main function
if __name__ == "__main__":
    main()
