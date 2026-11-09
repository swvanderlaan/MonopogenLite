#!/bin/bash

#SBATCH --job-name=mpg_germline       # Job name
#SBATCH --output=/sfs/gpfs/tardis/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen/mpg_germline_%A_%a.out  # Standard output log
#SBATCH --error=/sfs/gpfs/tardis/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen/mpg_germline_%A_%a.err   # Error log
#SBATCH --array=1-23   # One task per chromosome: 1-22, and 23 = X
#SBATCH --ntasks=1     # Number of tasks
#SBATCH --cpus-per-task=8  # Number of CPU cores per task; matches --nthreads below
#SBATCH --mem=32G                    # Memory per node (specify in GB)
#SBATCH --time=12:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl      # Where to send mail
#SBATCH --partition=standard # specific to rivanna
#SBATCH -A cphg-millerlab # specific to rivanna

# Description:
# Run germline variant calling of MonopogenLite for one study, one chromosome per SLURM array task.
# All samples of the study are called together per chromosome. The per-chromosome BAM lists are
# built by MonopogenLite from the preProcess output (one monopogen_<sample> folder per sample).

# --- Define script metadata ---
VERSION="v1.3.0"
VERSION_DATE="2026-11-09"
VERSION_NAME="runGermline"
VERSION_NAME_TEXT="Calling SNPs for scRNA-seq or snATAC-seq data using MonopogenLite."
COPYRIGHT="Copyright 1979-2026. José Verdezoto Mosquera; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science."
COPYRIGHT_TEXT="This script is licensed under the MIT license."

# --- Input arguments ---
STUDY=""
DRY_RUN=false

# --- Define default paths (same defaults as runPreprocess.sh) ---
MPG="/sfs/gpfs/tardis/project/cphg-millerlab/software/MonopogenLite"
PROJECT_DIR="/sfs/gpfs/tardis/project/cphg-millerlab/MetaPlaq/v2/processed_data/sc_germline/monopogen"
# MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"
# PROJECT_DIR="/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq"

# --- Reference genome, GRCh38 ---
GRCh38="/project/cphg-millerlab/Jose/Genome_assets/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

display_help() {
    echo ""
    echo "========================================================================"
    echo "$VERSION_NAME $VERSION ($VERSION_DATE)"
    echo "$VERSION_NAME_TEXT."
    echo "========================================================================"
    echo ""
    echo "Usage:"
    echo "sbatch runGermline.sh --study <name> [--input-dir <dir>] [--dry-run|-d] [--help|-h]"
    echo ""
    echo "Description:"
    echo "This script runs the germline variant calling of MonopogenLite for one study on a"
    echo "SLURM HPC cluster (e.g. Rivanna). Each array task (1-23) calls one chromosome"
    echo "(23 = X) for all samples of the study together."
    echo ""
    echo "Arguments:"
    echo "  --study NAME            Name of the study (folder under --input-dir holding the monopogen_<sample> folders). Required."
    echo "  --input-dir DIR         Directory holding the study folders, i.e. the --project-dir of runPreprocess.sh. Default: \$PROJECT_DIR"
    echo "  --mpg DIR               Path to the MonopogenLite root directory. Default: $MPG"
    echo "  --project-dir DIR       Project directory. Default: $PROJECT_DIR"
    echo "  --reference FILE        Reference genome FASTA file. Default: $GRCh38"
    echo "  --imputation-panel DIR  Imputation panel directory. Default: \$MPG/resources/"
    echo "  --region-file FILE      Path to the region file. Default: \$PROJECT_DIR/region.lst"
    echo "  --dry-run, -d           Only write the job script for the chromosome, without running it."
    echo "  --help, -h              Show this help message and exit."
    echo ""
    echo "Environment:"
    echo "Requires conda or mamba environment named 'monopogen' with all dependencies loaded."
    echo ""
    echo "Notes:"
    echo "- This script is meant to be submitted via sbatch; the array is fixed at 1-23."
    echo "- Run runPreprocess.sh for all samples of the study first."
    echo ""
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
    echo "========================================================================"
}

# --- Parse command line arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --help|-h) display_help; exit 0;;
        --study)
            STUDY="$2"
            shift
            ;;
        --input-dir)
            INPUT_DIR="$2"
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

# --- Defaults that depend on other arguments ---
INPUT_DIR="${INPUT_DIR:-$PROJECT_DIR}"
IMP_PANEL="${IMP_PANEL:-${MPG}/resources/}"
REGION_FILE="${REGION_FILE:-${PROJECT_DIR}/region.lst}"
STUDY_DIR="${INPUT_DIR%/}/${STUDY}"

echo "============================================"
echo "Monopogen: Running germline variant calling"
echo "============================================"
echo ""

# --- Validate the input arguments ---
echo "Checking guardians for the provided arguments."
if [[ -z "$STUDY" ]]; then
    echo "Error: --study is required." >&2
    exit 1
fi
if [[ ! -d "$STUDY_DIR" ]]; then
    echo "Error: Study directory not found at ${STUDY_DIR}" >&2
    exit 1
fi
if [[ ! -f "$REGION_FILE" ]]; then
    echo "Error: Region file not found at ${REGION_FILE}" >&2
    exit 1
fi
if [[ -z "$SLURM_ARRAY_TASK_ID" || "$SLURM_ARRAY_TASK_ID" -lt 1 || "$SLURM_ARRAY_TASK_ID" -gt 23 ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID:-unset}) must be 1-23; submit with sbatch --array=1-23." >&2
    exit 1
fi

# --- Chromosome for this array task ---
if [[ "$SLURM_ARRAY_TASK_ID" -eq 23 ]]; then
    CHROMOSOME="X"
else
    CHROMOSOME="$SLURM_ARRAY_TASK_ID"
fi

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

echo "Checking variables:"
echo "> Monopogen path.................: ${MPG}"
echo "> bcftools version...............: $(bcftools --version | head -1)"
echo "> java version...................: $(java -version 2>&1 | head -1)"
echo ""
echo "> Study directory................: ${STUDY_DIR}"
echo "> Study selected.................: ${STUDY}"
echo "> Chromosome.....................: ${CHROMOSOME}"
echo "> Reference genome...............: ${GRCh38}"
echo "> Imputation panel...............: ${IMP_PANEL}"
echo "> Region file....................: ${REGION_FILE}"
echo "> Samples (preProcess folders)...: $(find "$STUDY_DIR" -mindepth 1 -maxdepth 1 -type d -name 'monopogen_*' | wc -l)"
if [ "$DRY_RUN" = true ]; then
    echo "> Dry run mode...................: enabled"
else
    echo "> Dry run mode...................: disabled"
fi

echo ""
echo "Running Monopogen"
echo "> Germline analysis for chromosome ${CHROMOSOME} in study ${STUDY}..."

MONOPOGEN_CMD=(python "${MPG}/src/MonopogenLite.py" germline
    --region "${REGION_FILE}"
    --reference "${GRCh38}"
    --imputation-panel "${IMP_PANEL}"
    --step all
    --app-path "${MPG}/apps"
    --nthreads 8
    --out "${STUDY_DIR}"
    --chromosome "${CHROMOSOME}"
    --verbose)
if [ "$DRY_RUN" = true ]; then
    MONOPOGEN_CMD+=(--norun)
fi

echo "Running command: [ ${MONOPOGEN_CMD[*]} ]"
"${MONOPOGEN_CMD[@]}"
STATUS=$?

# --- Deactivate the conda environment ---
if command -v micromamba &> /dev/null; then
    micromamba deactivate
else
    mamba deactivate
fi

if [ $STATUS -ne 0 ]; then
    echo "Error: germline variant calling failed for chromosome ${CHROMOSOME} in study ${STUDY}." >&2
    exit $STATUS
fi

echo "Wow. That was a lot. Let's have a beer, buddy!"
