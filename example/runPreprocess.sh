#!/bin/bash

#SBATCH --job-name=mpg_preprocess       # Job name
#SBATCH --output=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_preprocess_%A_%a.out     # Standard output log file
#SBATCH --error=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_preprocess_%A_%a.err   # Error log
#SBATCH --array=0-49   # Array range (adjust based on the size of SAMPLE_LIST)
#SBATCH --ntasks=1     # Number of tasks
#SBATCH --cpus-per-task=8  # Number of CPU cores per task
#SBATCH --mem=8G                    # Memory per node (specify in GB)
#SBATCH --time=01:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl      # Where to send mail

# Description: 
# Run preProcess of Monopogen on the provided BAM files.
# It will try to run all the samples sequentially using the --mem and --time 
# provided in the SLURM header.
# 
# Change log:
# v1.1.1. 2024-09-18. Changed reference to MonopogenLite.
# v1.0.0. 2024-09-11. Initial version.
# Version: 1.1.1
# Author: Sander W. van der Laan
# Date: 2024-09-18
# Usage: sbatch runPreprocess.sh
# Arguments:
#   --help, -h  Show this help message

# Setting some variables for the script
# macOS
# MPG="$HOME/git/Monopogen"
# HPC
MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"
PROJECT_DIR="/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq"

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

# Argument validation and help
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: sbatch runPreprocess.sh" >&2
    exit 0
elif [ $# -gt 0 ]; then
    echo "Invalid argument. Usage: sbatch runPreprocess.sh" >&2
    exit 1
fi

echo "Monopogen: Preprocessing BAM files"

# Load the required conda environment and check if conda activate was successful
echo "Loading required mamba environment containing the monopogen installation..."
eval "$(conda shell.bash hook)"
conda activate monopogen

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate monopogen environment." >&2
    exit 1
fi
echo "> Checking existence of relevant apps..."
samtools --version
bcftools --version
vcftools --version
bgzip --version
java -version

echo ""
echo "> running MonopogenLite"

# Get the sample for the current array task
SAMPLE=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}
echo "  - preProcessing Monopogen for sample ${SAMPLE}..."

# Define the file name for each sample (e.g., SAMPLE.bam.lst)
FILENAME="${SAMPLE}.bam.lst"

# Create the line with sample-specific information
LINE="${SAMPLE},/hpc/dhl_ec/data/_ae_originals/AESCRNA/raw/post_alignment_output/${SAMPLE}/small_sam/item-1/${SAMPLE}_merged.bam"

# Write the line to the file
echo "$LINE" > "${PROJECT_DIR}/monopogen/${SAMPLE}.bam.lst"

# Print confirmation
echo "  - file $FILENAME has been created..."

# Run MonopogenLite preProcess for this sample
echo "  - running MonopogenLite preProcess for sample ${SAMPLE}..."
python ${MPG}/src/MonopogenLite.py preProcess \
    --bamFile ${PROJECT_DIR}/monopogen/${SAMPLE}.bam.lst --out ${PROJECT_DIR}/monopogen/monopogen_${SAMPLE} \
    --app-path ${MPG}/apps \
    --max-mismatch 3 \
    --nthreads 8 \
    --platform-library celseq2 --verbose

echo "Wow. That was a lot. Let's have a beer, buddy!"

# Deactivate the conda environment
conda deactivate

# The MonopogenLite preProcess script has the following options:

# python src/MonopogenLite.py  preProcess --help
# usage: MonopogenLite.py preProcess [-h] -b BAMFILE [-o OUT] -a APP_PATH [-m MAX_MISMATCH] [-t NTHREADS] -l {10x,smartseq2,celseq2}
#                                    [-v] [-d]

# Preprocess of BAM files including removing reads with high alignment mismatches. Default mismatch threshold is 3.

# optional arguments:
#   -h, --help            show this help message and exit
#   -b BAMFILE, --bamFile BAMFILE
#                         The comma-separated listf of bam-files for the study sample. The first column should have the sampleID, and
#                         the second column the location of the corresponding bam-file. The bam-files should be sorted and indexed.
#                         If there are multiple samples, each row with each sample. Required. (default: None)
#   -o OUT, --out OUT     The output directory. The output will be saved in the output directory. Required. (default: None)
#   -a APP_PATH, --app-path APP_PATH
#                         The app library paths used in the tool. The app library paths should include (a symlink to) beagle. Also
#                         see wiki for installation instructions of relevant tools (samtools, bcftools, bgzip, and java). Required.
#                         (default: None)
#   -m MAX_MISMATCH, --max-mismatch MAX_MISMATCH
#                         The maximal alignment mismatch allowed in one reads for variant calling. Default is 3. (default: 3)
#   -t NTHREADS, --nthreads NTHREADS
#                         Number of threads used for SNVs calling. Default is 1. (default: 1)
#   -l {10x,smartseq2,celseq2}, --platform-library {10x,smartseq2,celseq2}
#                         The platform library used for sequencing. This can be 10x, smartseq2, or celseq2. Required. (default: None)
#   -v, --verbose         Increase output verbosity. (default: False)
#   -d, --debug           For debugging, specifically for installed tools. (default: False)