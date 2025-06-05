#!/bin/bash

SCRIPT="runGermline.sh"

display_help() {
  echo ""
  echo "Usage: $0 --input-dir DIR"
  echo ""
  echo "This script submits germline variant calling jobs for multiple studies."
  echo ""
  echo "Arguments:"
  echo "  --input-dir DIR     Path to the input directory containing study folders"
  echo "  --studies LIST       Comma-separated list of study names to process (overrides default list)"
  echo "  --help, -h          Display this help message and exit"
  echo ""
}

INPUT_DIR=""
STUDIES=()

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

if [ -z "$INPUT_DIR" ]; then
  echo "Error: --input-dir is required."
  display_help
  exit 1
fi

# Default list of studies if --studies is not provided
if [ ${#STUDIES[@]} -eq 0 ]; then
  STUDIES=(
    YOUR_STUDY_1
    YOUR_STUDY_2
    YOUR_STUDY_3
  )
fi

for STUDY in "${STUDIES[@]}"; do
  echo "Submitting $STUDY >>>"
  # Step 1: Get the array range
  ARRAY_RANGE=$(bash $SCRIPT --input-dir "$INPUT_DIR" --study "$STUDY" --generate-array | awk '{print $5}')

  if [ -z "$ARRAY_RANGE" ]; then
    echo "Warning: no array range generated for $STUDY. Skipping."
    continue
  fi
  echo "Array range for $STUDY: $ARRAY_RANGE"
  
  # Step 2: Submit the job
  sbatch $ARRAY_RANGE $SCRIPT --input-dir "$INPUT_DIR" --study "$STUDY"
done