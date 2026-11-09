#!/bin/bash

#SBATCH --job-name=mpg_preprocess       # Job name
#SBATCH --output=/sfs/gpfs/tardis/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen/mpg_preprocess_%A_%a.out     # Standard output log file
#SBATCH --error=/sfs/gpfs/tardis/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen/mpg_preprocess_%A_%a.err   # Error log
#SBATCH --array=0-49   # Array range (adjust based on the size of SAMPLE_LIST)
#SBATCH --ntasks=1     # Number of tasks
#SBATCH --cpus-per-task=8  # Number of CPU cores per task
#SBATCH --mem=8G                    # Memory per node (specify in GB)
#SBATCH --time=08:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl      # Where to send mail
#SBATCH --partition=standard # specific to rivanna
#SBATCH -A cphg-millerlab # specific to rivanna

# Description: 
# Run preProcess of Monopogen on the provided BAM files.
# It will try to run all the samples sequentially using the --mem and --time 
# provided in the SLURM header.

# --- Define script metadata ---
VERSION="v1.2.8"
VERSION_DATE="2026-01-28"
VERSION_NAME="runPreprocess"
VERSION_NAME_TEXT="Preprocess scRNA-seq or snATAC-seq data for usage with MonopogenLite."
COPYRIGHT="Copyright 1979-2026. José Verdezoto Mosquera; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science."
COPYRIGHT_TEXT="This script is licensed under the MIT license."

# --- Define default paths ---
MPG="/sfs/gpfs/tardis/project/cphg-millerlab/software/MonopogenLite"
PROJECT_DIR="/sfs/gpfs/tardis/project/cphg-millerlab/MetaPlaq/v2/processed_data/sc_germline/monopogen"
# MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"
# PROJECT_DIR="/hpc/dhl_ec/svanderlaan/projects/molqtl_scrnaseq"

display_help() {
    echo ""
    echo "========================================================================"
    echo "$VERSION_NAME $VERSION ($VERSION_DATE)"
    echo "$VERSION_NAME_TEXT."
    echo "========================================================================"
    echo ""
    echo "Usage:"
    echo "sbatch runPreprocess.sh --input-dir <dir> --study <name> --generate-array [--dry-run|-d] [--help|-h]"
    echo ""
    echo "Description:"
    echo "This script runs the preProcess of MonopogenLite on the provided BAM files."
    echo "It will try to run all the samples sequentially using the --mem and --time provided in the SLURM header." 
    echo ""
    echo "Arguments:"
    echo "  --input-dir <dir>   Path to the directory containing the study samples. REQUIRED."
    echo "  --study <name>      Name of the study to be processed. REQUIRED."
    echo "  --platform <name>   Platform library used for sequencing. REQUIRED. Options: 10x, smartseq2, celseq2."
    echo "  --mpg <dir>         Path to the MonopogenLite root directory. Default: $MPG"
    echo "  --project-dir <dir> Project directory with region.lst. Default: $PROJECT_DIR"
    echo "  --generate-array    Only output the recommended --array setting based on number of samples, then exit."
    echo "  --dry-run, -d       Perform a dry run to show what would be executed without running it."
    echo "  --help, -h          Show this help message and exit."
    echo ""
    echo "Environment:"
    echo "Requires conda or mamba environment named 'monopogen' with all dependencies loaded."
    echo ""
    echo "Notes:"
    echo "- This script is meant to be submitted via sbatch."
    echo "- It automatically calculates the job array length based on the number of samples."
    echo ""
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
    echo "========================================================================"
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
        --platform)
            PLATFORM="$2"
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

# --- Validate input arguments ---
validate_arguments() {
if [ -z "$INPUT_DIR" ] || [ -z "$STUDY" ] || [ -z "$PLATFORM" ]; then
    echo "Error: Missing input arguments. --input-dir, --study, and --platform are required." >&2
    display_help
    exit 1
fi
}

echo "============================================"
echo "Monopogen: Running preprocess variant calling"
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

# --- Display progress ---
echo "Monopogen: Preprocessing BAM files"

# --- Display progress ---
echo ""
echo "Running MonopogenLite"

# --- Define output directory ---
OUTPUT_DIR="${PROJECT_DIR}/${STUDY}"
mkdir -p "$OUTPUT_DIR"

# --- Optionally calculate and output correct array range ---
if [[ "$GENERATE_ARRAY" = true ]]; then
    # --- Get list of sample directories ---
    SAMPLES=($(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d))
    if [[ ${#SAMPLES[@]} -eq 0 ]]; then
        echo "Error: No sample directories found for study $STUDY." >&2
        exit 1
    fi
    echo "Suggested sbatch array range: --array=0-$(( ${#SAMPLES[@]} - 1 ))"
    exit 0
# --- If not generating array, proceed with processing ---
else
    echo "Running MonopogenLite preProcess for study: $STUDY"
    # --- Get list of sample directories ---
    SAMPLES=($(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d))
    if [ "${#SAMPLES[@]}" -eq 0 ]; then
        echo "Error: No sample directories found for study $STUDY." >&2
        exit 1
    fi

fi

# --- Validate SLURM_ARRAY_TASK_ID ---
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#SAMPLES[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) out of range." >&2
    exit 1
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_NAME=$(basename "$SAMPLE")

echo "> searching for and listing BAM-files for sample: $SAMPLE_NAME"

# --- Locate BAM file ---
# This is a more robust way to find the BAM file
# This assumes the BAM file is named in a specific way. Adjust the pattern as needed.
# This will find the first BAM file that matches the pattern in the sample's outs directory.
# BAM_FILE=$(find "${SAMPLE}/outs" -maxdepth 1 -type f -name "*possorted_*_bam.bam" | head -n 1)
# This method is more robust and deterministic, and explicitly checks for a specific file name.

# This is for most studies; 10X Genomics
if [ -f "${SAMPLE}/outs/possorted_genome_bam.bam" ]; then
    BAM_FILE="${SAMPLE}/outs/possorted_genome_bam.bam"
    echo "BAM-file found: [${BAM_FILE}]"
# This is for Turner et al. 2022; Mosquera_Auguste_et_al; Mosquera_et_al; 10X Genomics
elif [ -f "${SAMPLE}/outs/possorted_bam.bam" ]; then
    BAM_FILE="${SAMPLE}/outs/possorted_bam.bam"
    echo "BAM-file found: [${BAM_FILE}]"
# This is for Vacante et al.; 10X Genomics
elif [ -f "${SAMPLE}/outs/gex_possorted_bam.bam" ]; then
    BAM_FILE="${SAMPLE}/outs/gex_possorted_bam.bam"
    echo "BAM-file found: [${BAM_FILE}]"
# This is for Vacante et al.; Turner_et_al_2022; Amrute_WashU_multiome; 10X Genomics
elif [ -f "${SAMPLE}/outs/atac_possorted_bam.bam" ]; then
    BAM_FILE="${SAMPLE}/outs/atac_possorted_bam.bam"
    echo "BAM-file found: [${BAM_FILE}]"
# This for Katyayani et al.; STARsolo
elif [ -f "${SAMPLE}/${SAMPLE_NAME}_Aligned_sortedByCoord_updated_header_out.bam" ]; then
    BAM_FILE="${SAMPLE}/${SAMPLE_NAME}_Aligned_sortedByCoord_updated_header_out.bam"
    echo "BAM-file found: [${BAM_FILE}]"
else
    echo "Error: No valid BAM file (should be [gex_]possorted_*bam.bam) found in ${SAMPLE}[/outs]" >&2
    exit 1
fi

# Check if BAM file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found at expected location: $BAM_FILE" >&2
    exit 1
fi

# --- Create .bam.lst ---
BAM_LIST_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.bam.lst"

# --- Create the line with sample-specific information ---
LINE="${SAMPLE_NAME},${BAM_FILE}"

# --- Write the line to the file ---
echo "$LINE" > "${BAM_LIST_FILE}"

# --- Run MonopogenLite preProcess ---
echo "> preprocessing sample: $SAMPLE_NAME"
LOG_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_preprocess.log"
python ${MPG}/src/MonopogenLite.py preProcess \
    --bamFile "$BAM_LIST_FILE" \
    --out "${OUTPUT_DIR}/monopogen_${SAMPLE_NAME}" \
    --app-path "${MPG}/apps" \
    --max-mismatch 3 \
    --nthreads 8 \
    --platform-library "$PLATFORM" --verbose &> "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "Error processing sample ${SAMPLE_NAME}. See log: [${LOG_FILE}]." >&2
    exit 1
fi

echo "Preprocessing complete for ${SAMPLE_NAME} from ${STUDY} (in ${INPUT_DIR}). Output: [${OUTPUT_DIR}/monopogen_${SAMPLE_NAME}]."

echo "Wow. That was a lot. Let's have a beer, buddy!"

# --- Deactivate the conda environment ---
if command -v micromamba &> /dev/null; then
    micromamba deactivate
else
    mamba deactivate
fi

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
