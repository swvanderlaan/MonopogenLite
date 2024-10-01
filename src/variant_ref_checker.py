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
    command = ['plot-vcfstats', '--prefix', os.path.join(output_dir, 'plot'), '--main-title', main_title, stats_file]
    
    logger.debug(f"> Running `bcftools plot-vcfstats` with output directory [{output_dir}].")
    
    # Run the command
    # subprocess.run(command, check=True)
    result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logger.error(f"Error running [{command}]: {result.stderr.decode('utf-8')}.")
        raise RuntimeError(f"Oh oh, `bcftools plot-vcfstats` failed for [{input_vcf}].")

    # Move the stats file to the output directory
    logger.debug(f"> Moving stats-file to output directory [{output_dir}].")
    new_stats_file = os.path.join(output_dir, os.path.basename(stats_file))
    shutil.move(stats_file, new_stats_file)
    
    logger.debug(f"> Moved stats file to [{new_stats_file}].")
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

# Plot the allele frequency histogram using cmcrameri colormap
def plot_allele_frequency_histogram(allele_frequencies, snp_counts, output_prefix, chromosome, logger, verbose=False):
    """Plot the allele frequency distribution as a histogram using cmcrameri colormap."""
    if not allele_frequencies or not snp_counts:
        logger.warning("No data available for plotting allele frequency histogram.")
        return
    
    logger.debug(f"> Plotting allele frequency histogram for [{len(allele_frequencies)}] entries.")
    
    plt.figure(figsize=(10, 6))
    
    # Use cmcrameri colormap for the bar colors
    cmap = cmc.lipari
    norm = plt.Normalize(vmin=min(allele_frequencies), vmax=max(allele_frequencies))
    colors = cmap(norm(allele_frequencies))

    plt.bar(allele_frequencies, snp_counts, width=0.0005, color=colors, edgecolor='black')

    plt.title(f'Allele frequency chromosome {chromosome}')
    plt.xlabel(f'allele frequency')
    plt.ylabel(f'number of variants')

    # Add colorbar to show the color mapping
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, label='Allele Frequency')
    
    # Save the plot
    logger.debug(f"> Saving allele frequency histogram as PNG and PDF.")
    histogram_file = f"{output_prefix}_allele_frequency_histogram.png"
    plt.savefig(histogram_file)
    histogram_file_pdf = f"{output_prefix}_allele_frequency_histogram.pdf"
    plt.savefig(histogram_file_pdf)
    plt.close()
    
    logger.debug(f"> Allele frequency histogram saved as [{histogram_file} and {histogram_file_pdf}].")

# Generate custom plots
def generate_plots(input_vcf, output_dir, logger, verbose=False):
    """Generate custom plots using cmcrameri color maps."""

    logger.debug(f"> Generating custom plots with cmcrameri color maps in [{output_dir}].")

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
    logger.debug(f"> Saving plot as PNG and PDF in [{output_dir}].")
    plt.savefig(os.path.join(output_dir, f"{base_name}_sample_plot.png"))
    plt.savefig(os.path.join(output_dir, f"{base_name}_sample_plot.pdf"))

    logger.debug(f"> Plots saved as PNG and PDF in {output_dir}.")

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
    
    input_vcf = args.input
    chromosome = args.chr
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
    logger.info(f"{'='*100}\n")
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE})\n")
    logger.info(f"Settings")
    logger.info(f"Input file..........: {input_vcf}")
    logger.info(f"Chromosome..........: {chromosome}")
    logger.info(f"Output directory....: {output_dir}")
    logger.info(f"Output stats file...: {output_stats}")
    logger.info(f"Verbose.............: {args.verbose}\n")
    
    # Run bcftools stats
    logger.info(f"Calculating statistics for chromosome {chromosome}...")
    run_bcftools_stats(input_vcf, output_stats, logger, args.verbose)
    
    # Run bcftools plot-vcfstats
    logger.info(f"Generating plots for chromosome {chromosome}...")
    run_bcftools_plot(output_stats, output_prefix, chromosome, logger, args.verbose)

    # Parse the stats file for allele frequency data
    logger.info(f"Parsing allele frequency data from {output_stats}...")
    allele_frequencies, snp_counts = parse_allele_frequency(output_stats, logger, args.verbose)
    
    # Plot the allele frequency distribution
    logger.info("Plotting allele frequency histogram...")
    plot_allele_frequency_histogram(allele_frequencies, snp_counts, output_prefix, chromosome, logger, args.verbose)
    
    # Generate additional custom plots
    # generate_plots(input_vcf, output_dir, args.verbose)

    logger.info(f"Log file saved as {datetime.now().strftime('%Y%m%d')}.variant_ref_checker.log")
    logger.info(f"Output files saved in {output_dir}.\n")

    logger.info("All done. Let's have a beer, buddy!")
    logger.info(f"{'='*100}\n")
    logger.info(f"{VERSION_NAME} v{VERSION} ({VERSION_DATE}) | {COPYRIGHT}.\n")
    logger.info(f"{'='*100}\n")
# Run the main function
if __name__ == "__main__":
    main()
