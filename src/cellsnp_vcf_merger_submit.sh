#!/bin/bash

#SBATCH --job-name=submit_cellsnp_vcf_merger
#SBATCH --output=/sfs/gpfs/tardis/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen/cellsnp_vcf_merger_submit_%A_%a.out
#SBATCH --error=/sfs/gpfs/tardis/project/cphg-millerlab/swvanderlaan_man7zh/MetaPlaq/monopogen/cellsnp_vcf_merger_submit_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH -A cphg-millerlab

# Description: 
# Submit cellsnp_vcf_merger jobs to Rivanna.
# It will try to run all the samples sequentially using the --mem and --time 
# provided in the SLURM header.
# 
# Change log:
# v1.2.0. 2025-06-03.
# v1.1.0. 2025-06-03.
# v1.0.12. 2025-06-03.
# v1.0.11. 2025-06-03.
# v1.0.10. 2025-06-03.
# v1.0.9. 2025-06-04. No changes, just tracking version with cellsnp_vcf_merger_submit.sh.
# v1.0.8. 2025-06-03.
# v1.0.7. 2025-06-02. 
# v1.0.6. 2025-06-02. 
# v1.0.5. 2025-06-02. No changes, just tracking version with cellsnp_vcf_merger.py.
# v1.0.4. 2025-06-02. No changes, just tracking version with cellsnp_vcf_merger.py.
# v1.0.3. 2025-06-02.
# v1.0.2. 2025-06-02. 
# v1.0.1. 2025-06-02. 
# v1.0.0. 2025-04-01. Initial version.
# Version: 1.2.0
# Author: Sander W. van der Laan
# Date: 2025-06-03
# Usage: bash submit_cellsnp_vcf_merger.sh

# --- Parse arguments ---
INPUT_DIR=""
STUDY_NAME=""
MERGE_ALL=false
VERBOSE=false
DEBUG=false
CLEAN_UP=false
DRY_RUN=false
EXEC_MODE=""

print_help() {
    echo "Usage: bash submit_cellsnp_vcf_merger.sh [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --input-dir DIR            Input directory containing study VCFs"
    echo "  --type [slurm|local]       Execution mode"
    echo "  --study-name NAME          Study name to process (mutually exclusive with --merge-all-studies)"
    echo "  --merge-all-studies        Merge all studies (mutually exclusive with --study-name)"
    echo ""
    echo "Optional:"
    echo "  --clean-up                 Clean up temporary files after processing"
    echo "  --verbose                  Enable verbose output"
    echo "  --debug                    Enable debug output"
    echo "  --dry-run                  Show commands without executing"
    echo "  --help                     Show this help message and exit"
    echo ""
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input-dir) INPUT_DIR="$2"; shift ;;
        --study-name) 
            if [[ "$MERGE_ALL" == true ]]; then
                echo "Error: --study-name and --merge-all-studies are mutually exclusive." >&2
                exit 1
            fi
            STUDY_NAME="$2"
            shift
            ;;
        --merge-all-studies)
            if [[ -n "$STUDY_NAME" ]]; then
                echo "Error: --merge-all-studies and --study-name are mutually exclusive." >&2
                exit 1
            fi
            MERGE_ALL=true
            ;;
        --verbose) VERBOSE=true ;;
        --debug) DEBUG=true ;;
        --clean-up) CLEAN_UP=true ;;
        --dry-run) DRY_RUN=true ;;
        --type)
            case "$2" in
                slurm|local) EXEC_MODE="$2" ;;
                *) echo "Error: --type must be 'slurm' or 'local'." >&2; exit 1 ;;
            esac
            shift
            ;;
        --help) print_help; exit 0 ;;
        *) echo "Unknown argument: $1"; print_help; exit 1 ;;
    esac
    shift
done

# --- Validate required arguments ---
if [[ -z "$INPUT_DIR" || -z "$EXEC_MODE" ]]; then
    echo "Error: --input-dir and --type are required." >&2
    exit 1
fi
if [[ "$MERGE_ALL" == false && -z "$STUDY_NAME" ]]; then
    echo "Error: One of --study-name or --merge-all-studies must be specified." >&2
    exit 1
fi

# --- Load conda environment ---
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

# --- Determine script path ---
if [[ "$EXEC_MODE" == "local" ]]; then
    SCRIPT_PATH="./cellsnp_vcf_merger.py"
else
    SCRIPT_PATH="/sfs/gpfs/tardis/project/cphg-millerlab/MetaPlaq/cellsnp/cellsnp_vcf_merger.py"
fi

if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Error: cellsnp_vcf_merger.py not found at $SCRIPT_PATH" >&2
    exit 1
fi

# --- Build command ---
CMD=(python3 "$SCRIPT_PATH" --input-dir "$INPUT_DIR")
[[ "$MERGE_ALL" == true ]] && CMD+=(--merge-all-studies)
[[ -n "$STUDY_NAME" ]] && CMD+=(--study-name "$STUDY_NAME" --output-dir "$INPUT_DIR/$STUDY_NAME")
[[ "$VERBOSE" == true ]] && CMD+=(--verbose)
[[ "$DEBUG" == true ]] && CMD+=(--debug)
[[ "$CLEAN_UP" == true ]] && CMD+=(--clean-up)

# --- Run or Dry-run ---
echo "Running in $EXEC_MODE mode. Dry-run: $DRY_RUN"
echo "Command: ${CMD[*]}"
if [[ "$DRY_RUN" == false ]]; then
    "${CMD[@]}"
fi

# --- Deactivate conda environment ---
if command -v micromamba &> /dev/null; then
    micromamba deactivate
elif command -v mamba &> /dev/null; then
    mamba deactivate
else
    echo "No known conda environment manager found to deactivate."
fi