#!/bin/bash


# Description: 
# Submit preProcess of Monopogen on the provided BAM files.
# It will try to run all the samples sequentially using the --mem and --time 
# provided in the SLURM header.
# 
# Change log:
# Version: 1.2.6
# Author: Sander W. van der Laan
# Date: 2025-06-05
# Usage: bash runPreprocess_submit_jobs.sh

# Define base directory where studies are stored
STUDY_DIR="/sfs/gpfs/tardis/project/cphg-millerlab/MetaPlaq/v2/raw_data/post_alignment_outs/cellranger_outs"

# List of available studies (can be manually updated if needed)
STUDIES=("YOUR_STUDY_1" "YOUR_STUDY_2" "YOUR_STUDY_3")

# Loop over studies and submit array jobs
for STUDY in "${STUDIES[@]}"; do
    STUDY_PATH="${STUDY_DIR}/${STUDY}"

    if [ -d "$STUDY_PATH" ]; then
        echo "Submitting SLURM job for study: $STUDY"

        # Get sample count for SLURM array job
        SAMPLE_COUNT=$(find "$STUDY_PATH" -mindepth 1 -maxdepth 1 -type d | wc -l)

        if [ "$SAMPLE_COUNT" -eq 0 ]; then
            echo "Warning: No samples found for study $STUDY. Skipping..."
            continue
        fi

        # Submit SLURM job with correct array size
        sbatch --array=0-$(($SAMPLE_COUNT - 1)) runPreprocess_rivanna.sh "$STUDY_PATH" "$STUDY"
    else
        echo "Warning: Study directory $STUDY_PATH does not exist. Skipping..."
    fi
done