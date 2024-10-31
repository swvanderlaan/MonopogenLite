#!/bin/bash

#SBATCH --job-name=mpg_germline       # Job name
#SBATCH --output=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_germline_%A_%a.out  # Standard output and error log
#SBATCH --error=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_germline_%A_%a.err   # Error log
#SBATCH --array=0-49   # Array range (adjust based on the size of SAMPLE_LIST)
#SBATCH --ntasks=1     # Number of tasks
#SBATCH --cpus-per-task=8  # Number of CPU cores per task
#SBATCH --mem=32G                    # Memory per node (specify in GB)
#SBATCH --time=03:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl      # Where to send mail

# Description: 
# Run germline variant calling of Monopogen on the provided BAM files.
# It will try to run all the samples sequentially using the --mem and --time
# provided in the SLURM header.
# 
# Change log:
# v1.1.2, 2024-09-19. Added optional --dry-run argument.
# v1.1.1, 2024-09-18. Changed reference to MonopogenLite.
# v1.0.0, 2024-09-11. Initial version.
# Version: 1.1.2
# Author: Sander W. van der Laan
# Date: 2024-09-19
# Usage: sbatch runGermline.sh
# Arguments:
#   --help, -h       Show this help message
#   --dry-run, -d    Do a dry run without executing the jobs

# Argument validation and help
DRY_RUN=false
for arg in "$@"; do
    case $arg in
        --help|-h)
            echo "Usage: sbatch runGermline.sh [--dry-run|-d]" >&2
            exit 0
            ;;
        --dry-run|-d)
            DRY_RUN=true
            ;;
        *)
            echo "Invalid argument: $arg. Use --help or -h for usage." >&2
            exit 1
            ;;
    esac
done

echo "============================================"
echo "Monopogen: Running germline variant calling"
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
# samtools --version
# bcftools --version
# vcftools --version
# bgzip --version
# java -version
# refgenie --version

echo "> setting some variables"
# macOS
# MPG="$HOME/git/Monopogen"
# HPC
MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"
PROJECT_DIR="/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq"

# reference genome, GRCh38
GRCh38="/hpc/dhl_ec/data/references/fasta/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

# imputation panel for the 1000G phase 3, high coverage, 2504 samples
IMP_PANEL="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/"

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
echo "> bcftools version...............: $(bcftools --version)"
echo "> refgenie version...............: $(refgenie --version)"
echo "> java version...................: $(java --version)"
echo ""
echo "> Project directory..............: ${PROJECT_DIR}"
echo "> Reference genome...............: ${GRCh38}"
echo "> Imputation panel...............: ${IMP_PANEL}"
echo "> Number of samples to process...: ${#SAMPLE_LIST[@]}"
if [ "$DRY_RUN" = true ]; then
    echo "> Dry run mode...................: enabled"
else 
    echo "> Dry run mode...................: disabled"
fi

echo ""
echo "Running Monopogen"
# Get the sample for the current array task
SAMPLE=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}

# Print confirmation
echo "> Germline analysis for sample ${SAMPLE}..."

# Construct the command based on whether it's a dry run or actual execution
if [ "$DRY_RUN" = true ]; then
    echo "Dry run mode activated. No jobs will be executed."
    MONOPOGEN_CMD="python ${MPG}/src/MonopogenLite.py germline \
        --region ${PROJECT_DIR}/monopogen/region.lst \
        --reference ${GRCh38} \
        --imputation-panel ${IMP_PANEL} \
        --step all \
        --max-softClipped 3 \
        --app-path ${MPG}/apps \
        --nthreads 8 \
        --out ${PROJECT_DIR}/monopogen/monopogen_${SAMPLE} --verbose --norun"
else
    MONOPOGEN_CMD="python ${MPG}/src/MonopogenLite.py germline \
        --region ${PROJECT_DIR}/monopogen/region.lst \
        --reference ${GRCh38} \
        --imputation-panel ${IMP_PANEL} \
        --step all \
        --max-softClipped 3 \
        --app-path ${MPG}/apps \
        --nthreads 8 \
        --out ${PROJECT_DIR}/monopogen/monopogen_${SAMPLE} --verbose"
fi

# Run the constructed command
echo "Running command: [ ${MONOPOGEN_CMD} ]"
$MONOPOGEN_CMD

echo "Wow. That was a lot. Let's have a beer, buddy!"

# Deactivate the conda environment
mamba deactivate

# The Monopogen germline variant calling script has the following options:

# python src/MonopogenLite.py  germline --help
# usage: MonopogenLite.py germline [-h] -r REGION -s {varScan,varProb,varPhasing,all} [-o OUT] -g REFERENCE -p
#                                  IMPUTATION_PANEL [-m MAX_SOFTCLIPPED] -a APP_PATH [-t NTHREADS] [-n] [-v] [-d]

# Perform germline variant calling and phasing from single-cell sequencing data.

# optional arguments:
#   -h, --help            show this help message and exit
#   -r REGION, --region REGION
#                         The genome regions for variant calling. This file should have either 1 column (chromosome) or 3
#                         columns (chromosome, start, end), where chromosome X is noted as chrX. Required. (default: None)
#   -s {varScan,varProb,varPhasing,all}, --step {varScan,varProb,varPhasing,all}
#                         Run germline variant calling step by step. varScan: variant calling; varProb: variant phased
#                         genotype probabilities; varPhasing: variant phasing; all: all steps. Default is all. (default:
#                         all)
#   -o OUT, --out OUT     The output directory. The output will be saved in the output directory. Required. (default: None)
#   -g REFERENCE, --reference REFERENCE
#                         The human genome reference used for alignment (default: None)
#   -p IMPUTATION_PANEL, --imputation-panel IMPUTATION_PANEL
#                         The population-level variant panel for variant phasing, such as 1000 Genome phase 3 high-coverage
#                         b38 data. (default: None)
#   -m MAX_SOFTCLIPPED, --max-softClipped MAX_SOFTCLIPPED
#                         The maximal soft-clipped allowed in one reads for variant calling (default: 1)
#   -a APP_PATH, --app-path APP_PATH
#                         The app library paths used in the tool (default: None)
#   -t NTHREADS, --nthreads NTHREADS
#                         Number of jobs used for SNVs calling (default: 1)
#   -n, --norun           Generate the job scripts only. The jobs will not be run. (default: False)
#   -v, --verbose         Increase output verbosity (default: False)
#   -d, --debug           For debugging, specifically for installed tools and some intermediate steps. (default: False)