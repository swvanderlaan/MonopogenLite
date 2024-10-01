#!/usr/bin/env python3

# Import required modules
import argparse
import os
import subprocess
import shutil
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc # needed for custom color maps generate_plots function
import logging
from datetime import datetime

# Change log:
# * v1.0.2 2024-10-01: Added plotting of allele frequency histogram. Expanded logging.
# * v1.0.1 2024-10-01: Added flow where a temporary directory is created for bcftools plot-vcfstats.
# * v1.0.0 2024-09-25: Initial version. 
# Version and license information 
VERSION_NAME = "Variant Reference Checker"
VERSION = "1.0.2"
VERSION_DATE = "2024-10-01"
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

# Define the UtrechtSciencePark color scheme as a dictionary
# Website to convert HEX to RGB: http://hex.colorrrs.com.
# For some functions you should divide these numbers by 255.
UTRECHT_COLOR_SCHEME = {
    1: "#FBB820",   # yellow
    2: "#F59D10",   # gold
    3: "#E55738",   # salmon
    4: "#DB003F",   # darkpink
    5: "#E35493",   # lightpink
    6: "#D5267B",   # pink
    7: "#CC0071",   # hardpink
    8: "#A8448A",   # lightpurple
    9: "#9A3480",   # purple
    10: "#8D5B9A",  # lavendel
    11: "#705296",  # bluepurple
    12: "#686AA9",  # purpleblue
    13: "#6173AD",  # lightpurpleblue
    14: "#4C81BF",  # seablue
    15: "#2F8BC9",  # skyblue
    16: "#1290D9",  # azurblue
    17: "#1396D8",  # lightazurblue
    18: "#15A6C1",  # greenblue
    19: "#5EB17F",  # seaweedgreen
    20: "#86B833",  # yellowgreen
    21: "#C5D220",  # lightmossgreen
    22: "#9FC228",  # mossgreen
    23: "#78B113",  # lightgreen (for X chromosome)
    24: "#49A01D",  # green (for Y chromosome)
    25: "#595A5C",  # grey (for XY)
    26: "#A2A3A4",  # lightgrey (for MT)
    # Additional colors if needed
    27: "#D7D8D7",  # midgrey
    28: "#ECECEC",  # very lightgrey
    29: "#FFFFFF",  # white
    30: "#000000"   # black
}

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

# Check if bcftools is installed
def check_tool_availability(tool_name, logger):
    """Check if a required tool is available on the system."""
    logger.debug(f"> Checking if tool '{tool_name}' is available.")
    if shutil.which(tool_name) is None:
        logger.error(f"> Required tool '{tool_name}' is not available on the system.")
        raise RuntimeError(f"Oh oh, required tool '{tool_name}' is missing.")
    else:
        logger.debug(f">> Tool '{tool_name}' is available.")

# Run bcftools stats -- cannot run shell-directions '>', '<', etc.
def run_bcftools_stats(input_vcf, output_stats, logger, verbose=False):
    """Run `bcftools stats` to generate statistics from the VCF file."""
    command = ['bcftools', 'stats', input_vcf]
    
    logger.debug(f"> Calculating statistics on [{input_vcf}].")

    with open(output_stats, 'w') as outfile:
        result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logger.error(f"Error running [{command}]: {result.stderr.decode('utf-8')}.")
            raise RuntimeError(f"Oh oh, `bcftools stats` failed for [{input_vcf}].")

    logger.debug(f"> Summary statistics saved to [{output_stats}].")

# Run bcftools plot-vcfstats
def run_bcftools_plot(stats_file, output_prefix, chromosome, logger, verbose=False):
    """Run `bcftools plot-vcfstats` to generate plots from the stats file and move the stats file."""
    # Create a directory for this job using the output_prefix
    output_dir = f"{output_prefix}_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Generate the main title dynamically based on the chromosome
    logger.debug(f"> Generating main title for chromosome [{chromosome}].")
    main_title = f"Summary chr{chromosome}"

    # Run plot-vcfstats with the directory we just created as the output and pass the main title
    command = ['plot-vcfstats', '--prefix', output_dir, '--main-title', main_title, stats_file]
    
    logger.debug(f"> Running `bcftools plot-vcfstats` with output directory [{output_dir}].")
    
    # Run the command
    # subprocess.run(command, check=True)
    result = subprocess.run(command, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logger.error(f"Error running [{command}]: {result.stderr.decode('utf-8')}.")
        raise RuntimeError(f"Oh oh, `bcftools plot-vcfstats` failed for [{stats_file}].")

    # # Move the stats file to the output directory
    # logger.debug(f"> Moving stats-file to output directory [{output_dir}].")
    # new_stats_file = os.path.join(output_dir, os.path.basename(stats_file))
    # shutil.move(stats_file, new_stats_file)
    
    # logger.debug(f"> Moved stats file to [{new_stats_file}].")
    logger.debug(f"> Output files generated and saved in [{output_dir}].")

# Parse the allele frequency data
def parse_allele_frequency(stats_file, logger, verbose=False):
    """Parse the stats file to extract allele frequency data and number of SNPs."""
    allele_frequencies = []
    snp_counts = []
    
    logger.debug(f"> Parsing allele frequency data from [{stats_file}].")
    with open(stats_file, 'r') as f:
        for line in f:
            if line.startswith("AF"):  # Only process lines that start with AF
                fields = line.strip().split('\t')
                allele_frequencies.append(float(fields[2]))  # [3]allele frequency
                snp_counts.append(int(fields[3]))  # [4]number of SNPs
    
    if not allele_frequencies:
        logger.warning(f"No allele frequency data found in [{stats_file}].")
    else:
        logger.debug(f"> Extracted [{len(allele_frequencies)}] allele frequency entries from [{stats_file}].")
    
    logger.debug(f"> Extracted allele frequencies.")
    return allele_frequencies, snp_counts

def plot_allele_frequency_histogram(allele_frequencies, snp_counts, output_prefix, chromosome, logger, verbose=False):
    """Plot the allele frequency distribution as a histogram using UtrechtSciencePark color scheme."""
    if not allele_frequencies or not snp_counts:
        logger.warning("No data available for plotting allele frequency histogram.")
        return
    
    logger.debug(f"> Plotting allele frequency histogram for [{len(allele_frequencies)}] entries.")
    
    # Create a new figure
    plt.figure(figsize=(10, 6))

    # Map chromosome 'X' to 23 and 'Y' to 24, 'XY' to 25, and 'MT' to 26, and use the Utrecht color scheme
    if chromosome == 'X':
        chromosome_int = 23
    elif chromosome == 'Y':
        chromosome_int = 24
    elif chromosome == 'XY':
        chromosome_int = 25
    elif chromosome == 'MT':
        chromosome_int = 26
    else:
        chromosome_int = int(chromosome)  # For numeric chromosomes
    # Get the color for the given chromosome from the Utrecht color scheme
    chromosome_color = UTRECHT_COLOR_SCHEME.get(int(chromosome), "#000000")  # Default to black if not found
    
    # Create bar plot with the specified color for the chromosome
    plt.bar(allele_frequencies, snp_counts, width=0.0005, color=chromosome_color, edgecolor='black')

    # Title and labels
    plt.title(f'Allele frequency chromosome {chromosome}')
    plt.xlabel(f'allele frequency')
    plt.ylabel(f'number of variants')

    # Save the plot to the output_prefix_plots directory
    output_dir = f"{output_prefix}_plots"
    
    # Extract the base name of output_prefix (without any directory components)
    output_base_name = os.path.basename(output_prefix)

    logger.debug(f"> Saving allele frequency histogram as PNG and PDF in [{output_dir}].")
    histogram_file = os.path.join(output_dir, f"{output_base_name}.AF_histogram.png")
    plt.savefig(histogram_file)
    histogram_file_pdf = os.path.join(output_dir, f"{output_base_name}.AF_histogram.pdf")
    plt.savefig(histogram_file_pdf)
    # Close the plot
    plt.close()
    
    logger.debug(f"> Allele frequency histogram saved as [{histogram_file} and {histogram_file_pdf}].")

# Main function
def main():
    parser = argparse.ArgumentParser(description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Run bcftools stats and generate plots using cmcrameri colors.""",
		epilog=f"""
For a given input VCF file (`--input-vcf`) and the associated chromosome 
(`--chr`), this script will calculate for the given chromosome some
summary statistics. Next, it will generate summary plots including a histogram
of the allele frequencies. The output files will be saved in a directory
with the same name as the input VCF file.

Optionally, you can use the `--verbose` flag to enable verbose output.

Example:
    python variant_ref_checker.py --input input.chr1.vcf.gz --chr 1 (--verbose)

+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} +
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', required=True, help="Path to the input VCF file (vcf.gz format). Required.")
    parser.add_argument('-c', '--chr', help="The chromosome analyzed; needed for plotting. Required.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output. Optional.")
    parser.add_argument('-V', '--version', action='version', version=f'{VERSION_NAME} {VERSION} ({VERSION_DATE})')
    
    args = parser.parse_args()
    
    # Check if the --input and --chr are provided
    if not args.input or not args.chr:
        logger.error("Please provide both the input VCF file (`--input`) and the chromosome (`--chr`).")
        raise ValueError("Oh, oh. Please provide both the input VCF file and the chromosome.")

    input_vcf = args.input
    chromosome = args.chr
    output_dir = os.path.dirname(input_vcf)
    
    # Generate the stats file and plot prefixes
    base_name = os.path.basename(input_vcf).replace(".vcf.gz", "")
    output_stats = os.path.join(output_dir, f'{base_name}.vcf.stats.txt')
    output_prefix = os.path.join(output_dir, base_name)

    # Ensure the input file is a .vcf.gz file
    if not input_vcf.endswith(".vcf.gz"):
        raise ValueError("Oh oh. Input file must be a .vcf.gz file.")
    
    # Set up logger
    logger = setup_logger("variant_ref_checker", args.verbose)
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})\n")
    logger.info(f"Settings")
    logger.info(f"Input file..........: {input_vcf}")
    logger.info(f"Chromosome..........: {chromosome}")
    logger.info(f"Output directory....: {output_dir}")
    logger.info(f"Output stats file...: {output_stats}")
    logger.info(f"Verbose.............: {args.verbose}\n")
    
    # Check if required tools are available
    logger.info(f"> Checking availability of required tools.")
    check_tool_availability("bcftools", logger)
    check_tool_availability("plot-vcfstats", logger)
    
    # Run bcftools stats
    logger.info(f"Calculating statistics...")
    run_bcftools_stats(input_vcf, output_stats, logger, args.verbose)
    
    # Run bcftools plot-vcfstats
    logger.info(f"Generating summary plots...")
    run_bcftools_plot(output_stats, output_prefix, chromosome, logger, args.verbose)

    # Parse the stats file for allele frequency data
    logger.info(f"Parsing allele frequency data from {output_stats}...")
    allele_frequencies, snp_counts = parse_allele_frequency(output_stats, logger, args.verbose)
    
    # Plot the allele frequency distribution
    logger.info(f"Plotting allele frequency histogram...")
    plot_allele_frequency_histogram(allele_frequencies, snp_counts, output_prefix, chromosome, logger, args.verbose)

    logger.info(f"Log file saved as {datetime.now().strftime('%Y%m%d')}.variant_ref_checker.log")
    logger.info(f"Output files saved in {output_dir}.\n")

    logger.info("All done. Let's have a beer, buddy!\n")
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}) | {COPYRIGHT}.\n")
# Run the main function
if __name__ == "__main__":
    main()
