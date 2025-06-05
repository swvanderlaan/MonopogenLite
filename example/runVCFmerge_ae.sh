#!/bin/bash

#SBATCH --job-name=mpg_vcf_merge       # Job name
#SBATCH --output=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_vcf_merge_%A_%a.out  # Standard output and error log
#SBATCH --error=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_vcf_merge_%A_%a.err   # Error log
#SBATCH --array=0-49   # Array range (adjust based on the size of SAMPLE_LIST)
#SBATCH --ntasks=1     # Number of tasks
#SBATCH --cpus-per-task=8  # Number of CPU cores per task
#SBATCH --mem=16G                    # Memory per node (specify in GB)
#SBATCH --time=08:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl      # Where to send mail

# Description: 
# Run VCF merging for germline variant calling. 
# It will try to run all the samples sequentially using the --mem and --time
# provided in the SLURM header.
# 
# Change log:
# * v1.0.1, 2024-10-10: Fixed issue with VCF filenames and order of operations.
# * v1.0.0, 2024-09-13: Initial version.
# Version: 1.0.1
# Author: Sander W. van der Laan
# Date: 2024-10-10
# Usage: sbatch runVCFmerge.sh
# Arguments:
#   -v, --verbose   Verbose mode on, optional (default: off)
#   -h, --help      Show this help message

version_name="runVCFmerge"
version="1.0.1"
version_date="2024-10-10"
version_copyright="Copyright 1979-2024. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science"
version_copyright_text="
The MIT License (MIT).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
"

# Argument parsing and help message
verbose="off"
while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -v|--verbose) verbose="off"; shift ;; # Verbose mode on, optional (default: off)
        -h|--help) # Show help message 
            echo "Usage: sbatch runVCFmerge.sh [--verbose|-v] [--help|-h]"
            echo "Arguments:"
            echo "  -v, --verbose   Verbose mode on, optional (default: off)."
            echo "  -h, --help      Show this help message."
            exit 0;;
        *) 
            echo "Invalid option: $1"
            exit 1;;
    esac
    shift
done

echo "============================================"
echo "Monopogen: merging VCF files"
echo "============================================"
echo ""
# Load the required conda environment and check if conda activate was successful
echo "Loading required mamba environment containing the monopogen installation..."
source ~/.bashrc
mamba activate monopogen

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate monopogen environment." >&2
    exit 1
fi
echo "> Checking existence of relevant apps..."
bcftools --version
tabix --version
refgenie --version

echo "> setting some variables"
# macOS
# MPG="$HOME/git/Monopogen"
# HPC
MPG="/hpc/local/Rocky8/dhl_ec/software/Monopogen"
PROJECT_DIR="/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen"
# reference genome, GRCh38
GRCh38=$(refgenie seek hg38/fasta)

# Samples to process
SAMPLE_LIST=("4432_UMC-DE-037_HYJWFBGX9" "4440_UMC-DE-005_AH32W2BGX9" "4443_UMC-DE-017_H33GYBGX9" 
             "4443_UMC-DE-018_H33GYBGX9" "4443_UMC-DE-019_H33GYBGX9" "4447_UMC-DE-016_H33GYBGX9" 
             "4447_UMC-DE-020_H33GYBGX9" "4447_UMC-DE-021_AH32W2BGX9" "4448_UMC-DE-041_HYJWFBGX9" 
             "4450_UMC-DE-022_AHT3MNBGX7" "4450_UMC-DE-023_AHT3MNBGX7" "4450_UMC-DE-024_AHT3MNBGX7" 
             "4452_UMC-DE-025_AHT3MNBGX7" "4452_UMC-DE-026_AHT3MNBGX7" "4452_UMC-DE-027_AHT3K3BGX7" 
             "4453_UMC-DE-028_AHT3K3BGX7" "4453_UMC-DE-029_AHT3K3BGX7" "4453_UMC-DE-030_AHT3K3BGX7" 
             "4455_UMC-DE-034_AH73J7BGX9" "4458_UMC-DE-033_AH73J7BGX9" "4459_UMC-DE-032_AH73J7BGX9" 
             "4459_UMC-DE-s345_HT73MBGXH" "4470_UMC-DE-031_AH73J7BGX9" "4472_UMC-DE-036_HYJWFBGX9" 
             "4477_UMC-DE-040_HYJWFBGX9" "4477_UMC-DE-s326_HT73MBGXH" "4478_UMC-DE-038_HYJWFBGX9" 
             "4478_UMC-DE-s329_HT73MBGXH" "4480_UMC-DE-039_HYJWFBGX9" "4486_UMC-DE-042_HYJWFBGX9" 
             "4487_UMC-DE-035_HYJWFBGX9" "4487_UMC-DE-100_HG2WNBGXB" "4488_UMC-DE-043_HYJWFBGX9" 
             "4489_UMC-DE-102_HG2WNBGXB" "4491_UMC-DE-105_HG2WNBGXB" "4495_UMC-DE-108_HG2WNBGXB" 
             "4496_UMC-DE-111_HG2WNBGXB" "4500_UMC-DE-072_HG2WNBGXB" "4500_UMC-DE-s310_HNYV2BGXH" 
             "4501_UMC-DE-075_HG2WNBGXB" "4502_UMC-DE-068_HG2WNBGXB" "4513_UMC-DE-084_HG2WNBGXB" 
             "4580_UMC-DE-s112_HT73MBGXH" "4587_UMC-DE-s115_HT73MBGXH" "4601_UMC-DE-s124_HT73MBGXH" 
             "4602_UMC-DE-s118_HT73MBGXH" "4605_UMC-DE-s121_HT73MBGXH" "4653_UMC-DE-s406_HT73MBGXH" 
             "4675_UMC-DE-s400_HT73MBGXH" "4676_UMC-DE-s403_HT3KLBGXH")

echo "Checking variables:"
echo "> Monopogen path.................: ${MPG}"
echo "> Project directory..............: ${PROJECT_DIR}"
echo "> Number of samples to process...: ${#SAMPLE_LIST[@]}"
echo ""
echo "Merging VCF files for germline analysis..."
# Get the sample for the current array task
SAMPLE=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}

# Print confirmation
echo "> processing sample ${SAMPLE}..."

# /hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/monopogen_4432_UMC-DE-037_HYJWFBGX9/germline
# Create a list of VCF files (assumed to be named like chr1.vcf.gz, chr2.vcf.gz, etc.)
vcf_files_gl=""
vcf_files_gp=""
vcf_files_gphased=""
if [ $verbose == "on" ]; then
    echo "> creating a list of VCF files - we use genotype probabilities (gp.vcf.gz) and genotypes (germline.vcf) data..."
fi
for CHROM in $(seq 1 22) X; do
    echo "  - Processing chromosome ${CHROM}..."
    
    gl_vcf="${PROJECT_DIR}/monopogen_${SAMPLE}/germline/chr${CHROM}.gl.vcf.gz"
    gp_vcf="${PROJECT_DIR}/monopogen_${SAMPLE}/germline/chr${CHROM}.gp.vcf.gz"
    gphased_vcf="${PROJECT_DIR}/monopogen_${SAMPLE}/germline/chr${CHROM}.phased.vcf.gz"

    normalized_gl_vcf="${PROJECT_DIR}/monopogen_${SAMPLE}/germline/chr${CHROM}.gl.norm.vcf.gz"
    normalized_gp_vcf="${PROJECT_DIR}/monopogen_${SAMPLE}/germline/chr${CHROM}.gp.norm.vcf.gz"
    normalized_gphased_vcf="${PROJECT_DIR}/monopogen_${SAMPLE}/germline/chr${CHROM}.phased.norm.vcf.gz"

    # First bgzip the genotypes VCF file
    echo "  > bgzipping the genotypes VCF file..."
    bgzip $gl_vcf
    bgzip $gp_vcf 

    echo "  > bgzipping the phased VCF file..."
    bgzip $gphased_vcf 

    # Next tabix the genotypes VCF file
    echo "  > indexing the genotypes VCF file..."
    tabix -p vcf $gl_vcf
    tabix -p vcf $gp_vcf

    # Indexing the phased VCF file
    echo "  > indexing the phased VCF file..."
    tabix -p vcf $gphased_vcf

    # Normalize the VCF files against the reference genome
    echo "  > Normalizing VCFs..."
    bcftools norm -f $GRCh38 -Oz -o $normalized_gl_vcf $gl_vcf
    bcftools norm -f $GRCh38 -Oz -o $normalized_gp_vcf $gp_vcf
    bcftools norm -f $GRCh38 -Oz -o $normalized_gphased_vcf $gphased_vcf

    # Index the normalized VCFs
    echo "  > Indexing the normalized VCFs..."
    tabix -p vcf $normalized_gl_vcf
    tabix -p vcf $normalized_gp_vcf
    tabix -p vcf $normalized_gphased_vcf

    # Append normalized files to the list for concatenation
    vcf_files_gl+=" $normalized_gl_vcf"
    vcf_files_gp+=" $normalized_gp_vcf"
    vcf_files_gphased+=" $normalized_gphased_vcf"
done

# Concatenate VCF files
if [ $verbose == "on" ]; then
    echo "> concatenating VCF files..."
fi
bcftools concat $vcf_files_gl -Oz -o ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gl.merged.vcf.gz
bcftools concat $vcf_files_gt -Oz -o ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gp.merged.vcf.gz
bcftools concat $vcf_files -Oz -o ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.phased.merged.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error: Failed to concatenate VCF files for sample ${SAMPLE}." >&2
    exit 1
fi

# Sort the concatenated VCF
if [ $verbose == "on" ]; then
    echo "> sorting the merged VCF..."
fi
bcftools sort ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gl.merged.vcf.gz -Oz -o ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gl.merged.sorted.vcf.gz
bcftools sort ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gp.merged.vcf.gz -Oz -o ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gp.merged.sorted.vcf.gz
bcftools sort ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.phased.merged.vcf.gz -Oz -o ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.phased.merged.sorted.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error: Failed to sort the merged VCF for sample ${SAMPLE}." >&2
    exit 1
fi

# Index the sorted and concatenated VCF
if [ $verbose == "on" ]; then
    echo "> indexing the sorted VCF..."
fi
tabix -p vcf ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gl.merged.sorted.vcf.gz
tabix -p vcf ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.gp.merged.sorted.vcf.gz
tabix -p vcf ${PROJECT_DIR}/monopogen_${SAMPLE}/germline/${SAMPLE}.phased.merged.sorted.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error: Failed to index the sorted VCF for sample ${SAMPLE}." >&2
    exit 1
fi

echo "VCF file merged, sorted, and indexed successfully for sample ${SAMPLE}."
echo "Wow. That was a lot. Let's have a beer, buddy!"

# Deactivate the conda environment
mamba deactivate