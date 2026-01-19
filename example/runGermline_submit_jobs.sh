#!/bin/bash

SCRIPT="runGermline_rivanna.sh"

# --- Function to display help message ---
display_help() {
  echo "Submit Germline Variant Calling Jobs"
  echo "==================================="
  echo ""
  echo "Usage: $0 --input-dir DIR"
  echo ""
  echo "This script submits germline variant calling jobs for multiple studies."
  echo ""
  echo "Arguments:"
  echo "  --input-dir DIR     Path to the input directory containing study folders"
  echo "  --studies LIST      Comma-separated list of study names to process (overrides default list)"
  echo "  --help, -h          Display this help message and exit"
  echo ""
  echo "Default studies if --studies is not provided:"
  echo "  Alsaigh_et_al_2022, Bashore_et_al_2024, Cheng_et_al, Chou_et_al_2021,"
  echo "  Chou_et_al_2022, Eberhardt_et_al_2023, Fernandez_et_al_2019, Jaiswal_et_al,"
  echo "  Katyayani_et_al, Paloschi_et_al, Pan_et_al_2020, Qian_et_al, Turner_et_al_2022,"
  echo "  Vacante_et_al, Wirka_et_al_2019"
  echo ""
  echo "Example:"
  echo "  $0 --input-dir /path/to/input --studies Alsaigh_et_al_2022,Bashore_et_al_2024"
  echo ""
  exit 0
}

# --- Main script starts here ---
INPUT_DIR=""
STUDIES=()

# --- Parse command line arguments ---
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --input-dir)
      INPUT_DIR="$2"
      shift
      ;;
    --studies)
      IFS=',' read -ra STUDIES <<< "$2"
      shift
      ;;
    --help|-h)
      display_help
      exit 0
      ;;
    *)
      echo "Unknown parameter passed: $1"
      display_help
      exit 1
      ;;
  esac
  shift
done

# --- Check if input directory is provided ---
if [ -z "$INPUT_DIR" ]; then
  echo "Error: --input-dir is required."
  display_help
  exit 1
fi

# --- Default list of studies if --studies is not provided ---
if [ ${#STUDIES[@]} -eq 0 ]; then
  STUDIES=(
    Alsaigh_et_al_2022
    Bashore_et_al_2024
    Cheng_et_al
    Chou_et_al_2021
    Chou_et_al_2022
    Eberhardt_et_al_2023
    Fernandez_et_al_2019
    Jaiswal_et_al
    Katyayani_et_al
    Paloschi_et_al
    Pan_et_al_2020
    Qian_et_al
    Turner_et_al_2022
    Vacante_et_al
    Wirka_et_al_2019
  )
fi

# --- Start loop to submit jobs for each study ---
for STUDY in "${STUDIES[@]}"; do
    echo "Submitting $STUDY germline variant calling job..."
    
    # Extract array range safely
    ARRAY_RANGE=$(bash "$SCRIPT" --input-dir "$INPUT_DIR" --study "$STUDY" --generate-array | awk '/--array/ {print $5}')

    if [ -z "$ARRAY_RANGE" ]; then
        echo "⚠️ Warning: No array range generated for $STUDY. Skipping."
        continue
    fi

    echo "Array range for $STUDY: $ARRAY_RANGE"
    sbatch "$ARRAY_RANGE" "$SCRIPT" --input-dir "$INPUT_DIR" --study "$STUDY"
done