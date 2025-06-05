#!/bin/bash

#SBATCH --job-name=mpg_germline       # Job name
#SBATCH --output=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_germline_%A_%a.out  # Standard output and error log
#SBATCH --error=/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq/monopogen/mpg_germline_%A_%a.err   # Error log
#SBATCH --array=0-0  # Placeholder; will be replaced dynamically if script is submitted with --generate-array
#SBATCH --ntasks=1     # Number of tasks
#SBATCH --cpus-per-task=8  # Number of CPU cores per task
#SBATCH --mem=32G                    # Memory per node (specify in GB)
#SBATCH --time=03:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl      # Where to send mail
#SBATCH --partition=standard
#SBATCH -A cphg-millerlab

# Description: 
# Run germline variant calling of Monopogen on the provided BAM files.

# --- Define script metadata ---
VERSION="v1.2.7"
VERSION_DATE="2025-06-05"
VERSION_NAME="runGermline"
VERSION_NAME_TEXT="Calling SNPs for scRNA-seq or snATAC-seq data using MonopogenLite."
COPYRIGHT="Copyright 1979-2025. José Verdezoto Mosquera; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science."
COPYRIGHT_TEXT="This script is licensed under the MIT license."

# --- Input arguments ---
INPUT_DIR=""
STUDY=""
DRY_RUN=false

# --- Define paths ---
MPG="/project/cphg-millerlab/software/MonopogenLite"
PROJECT_DIR="/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen"

# --- Reference genome, GRCh38 ---
GRCh38="/project/cphg-millerlab/Jose/Genome_assets/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

# --- Imputation panel for the 1000G phase 3, high coverage, 2504 samples ---
IMP_PANEL="${MPG}/resources/"

# --- Region file for the analysis ---
REGION_LIST="${PROJECT_DIR}/region.lst"

display_help() {
    echo ""
    echo "========================================================================"
    echo "$VERSION_NAME $VERSION ($VERSION_DATE)"
    echo "$VERSION_NAME_TEXT."
    echo "========================================================================"
    echo ""
    echo "Usage:"
    echo "sbatch runGermline_rivanna.sh --input-dir <dir> --study <name> --generate-array [--dry-run|-d] [--help|-h]"
    echo ""
    echo "Description:"
    echo "This script runs the germline variant calling pipeline from MonopogenLite"
    echo "for a given study on a SLURM HPC cluster (e.g. Rivanna). It loops through"
    echo ".bam.lst files and executes the full germline workflow for each sample."
    echo ""
    echo "Arguments:"
    echo "  --input-dir DIR         Path to the input directory containing sample folders with BAM lists."
    echo "  --study NAME            Name of the study to process (matches folder name under input-dir)."
    echo "  --mpg DIR               Path to the MonopogenLite root directory. Default: $MPG"
    echo "  --project-dir DIR       Project directory with region.lst. Default: $PROJECT_DIR"
    echo "  --reference FILE        Reference genome FASTA file. Default: $GRCh38"
    echo "  --imputation-panel DIR  Imputation panel directory. Default: $IMP_PANEL"
    echo "  --region-file FILE      Path to the region file. Default: $REGION_LIST"
    echo "  --generate-array        Only output the recommended --array setting based on number of samples, then exit."
    echo "  --dry-run, -d           Perform a dry run to show what would be executed without running it."
    echo "  --help, -h              Show this help message and exit."
    echo ""
    echo "Environment:"
    echo "Requires conda or mamba environment named 'monopogen' with all dependencies loaded."
    echo ""
    echo "Notes:"
    echo "- This script is meant to be submitted via sbatch."
    echo "- It automatically calculates the job array length based on the number of BAM list files."
    echo ""
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
    echo "========================================================================"
}

# --- Validate arguments ---
validate_arguments() {
# Check if required arguments are provided
if [[ -z "$INPUT_DIR" || -z "$STUDY" ]]; then
    echo "Error: --input-dir and --study are required." >&2
    exit 1
fi

# Check for region file using REGION_FILE if set, else REGION_LIST
if [ ! -f "${REGION_FILE:-$REGION_LIST}" ]; then
    echo "Error: Region file not found at ${REGION_FILE:-$REGION_LIST}" >&2
    exit 1
fi

}

# --- Parse command line arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --help|-h) display_help; exit 0;;
        --input-dir)
            INPUT_DIR="$2"
            shift
            ;;
        --study)
            STUDY="$2"
            shift
            ;;
        --mpg)
            MPG="$2"
            shift
            ;;
        --project-dir)
            PROJECT_DIR="$2"
            shift
            ;;
        --reference)
            GRCh38="$2"
            shift
            ;;
        --imputation-panel)
            IMP_PANEL="$2"
            shift
            ;;
        --region-file)
            REGION_FILE="$2"
            shift
            ;;
        --generate-array)
            GENERATE_ARRAY=true
            ;;
        --dry-run|-d)
            DRY_RUN=true
            ;;
        *)
            echo "Invalid argument: $1. Use --help or -h for usage." >&2
            exit 1
            ;;
    esac
    shift
done

# --- Optionally calculate and output correct array range ---
if [[ "$GENERATE_ARRAY" = true ]]; then
    MONOPOGEN_ROOT="${INPUT_DIR%/}/${STUDY}"
    SAMPLE_COUNT=$(find "$MONOPOGEN_ROOT" -maxdepth 1 -name "*.bam.lst" | wc -l)
    if [[ "$SAMPLE_COUNT" -eq 0 ]]; then
        echo "No .bam.lst files found in $MONOPOGEN_ROOT. Cannot generate array range." >&2
        exit 1
    fi
    echo "Suggested sbatch array range: --array=0-$((SAMPLE_COUNT - 1))"
    exit 0
fi

echo "============================================"
echo "Monopogen: Running germline variant calling"
echo "============================================"
echo ""
# --- Validate the input arguments ---
echo "Checking guardians for the provided arguments."
validate_arguments

# --- Load conda environment ---
echo "Loading necessary environment..."
source ~/.bashrc
if command -v micromamba &> /dev/null; then
    micromamba activate monopogen
elif command -v mamba &> /dev/null; then
    mamba activate monopogen
else
    echo "Error: Neither micromamba nor mamba is available in PATH." >&2
    exit 1
fi

# --- Check necessary tools ---
for tool in samtools bcftools vcftools bgzip java; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool is not installed or not in PATH." >&2
        exit 1
    fi
done

echo "> setting some variables"
# --- Samples to process ---
MONOPOGEN_ROOT="${INPUT_DIR%/}/${STUDY}"
SAMPLE_LIST=()
while IFS= read -r bamlist_file; do
    sample=$(basename "$bamlist_file" .bam.lst)
    SAMPLE_LIST+=("$sample")
done < <(find "$MONOPOGEN_ROOT" -maxdepth 1 -name "*.bam.lst")

ARRAY_LIMIT=$((${#SAMPLE_LIST[@]} - 1))

echo "Checking variables:"
echo "> Monopogen path.................: ${MPG}"
echo "> bcftools version...............: $(bcftools --version)"
echo "> refgenie version...............: $(refgenie --version)"
echo "> java version...................: $(java --version)"
echo ""
echo "> Monopogen root directory........: ${MONOPOGEN_ROOT}"
echo "> Project directory..............: ${PROJECT_DIR}"
echo "> Input directory................: ${INPUT_DIR}"
echo "> Study selected.................: ${STUDY}"
echo "> Reference genome...............: ${GRCh38}"
echo "> Imputation panel...............: ${IMP_PANEL}"
echo "> Region file....................: ${REGION_LIST}"
echo "> Number of samples to process...: ${#SAMPLE_LIST[@]}"
if [ "$DRY_RUN" = true ]; then
    echo "> Dry run mode...................: enabled"
else 
    echo "> Dry run mode...................: disabled"
fi

echo ""
echo "Running Monopogen"
# --- Get the sample for the current array task ---
SAMPLE=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}
if [[ "$SLURM_ARRAY_TASK_ID" -ge "${#SAMPLE_LIST[@]}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of samples (${#SAMPLE_LIST[@]})." >&2
    exit 1
fi

# --- Print confirmation ---
echo "> Germline analysis for sample ${SAMPLE} in study ${STUDY}..."

# --- Construct the command based on whether it's a dry run or actual execution ---
if [ "$DRY_RUN" = true ]; then
    echo "Dry run mode activated. No jobs will be executed."
    MONOPOGEN_CMD="python ${MPG}/src/MonopogenLite.py germline \
        --region ${REGION_FILE:-$REGION_LIST} \
        --reference ${GRCh38} \
        --imputation-panel ${IMP_PANEL} \
        --step all \
        --max-softClipped 3 \
        --app-path ${MPG}/apps \
        --nthreads 8 \
        --out ${MONOPOGEN_ROOT}/monopogen_${SAMPLE} --verbose --norun"
else
    MONOPOGEN_CMD="python ${MPG}/src/MonopogenLite.py germline \
        --region ${REGION_FILE:-$REGION_LIST} \
        --reference ${GRCh38} \
        --imputation-panel ${IMP_PANEL} \
        --step all \
        --max-softClipped 3 \
        --app-path ${MPG}/apps \
        --nthreads 8 \
        --out ${MONOPOGEN_ROOT}/monopogen_${SAMPLE} --verbose"
fi

# --- Run the constructed command ---
if [ "$DRY_RUN" = true ]; then
    echo "Dry run mode: would run the following command:"
    echo "${MONOPOGEN_CMD}"
else
    echo "Running command: [ ${MONOPOGEN_CMD} ]"
    $MONOPOGEN_CMD
fi

echo "Wow. That was a lot. Let's have a beer, buddy!"

# --- Deactivate the conda environment ---
if command -v micromamba &> /dev/null; then
    micromamba deactivate
else
    mamba deactivate
fi

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