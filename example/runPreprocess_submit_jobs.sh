#!/bin/bash

# Title: submit_jobs_runPreprocess_rivanna.sh
# Description: 
# Submit preProcess of Monopogen on the provided BAM files.
# It will try to run all the samples sequentially using the --mem and --time 
# provided in the SLURM header.
# 
# Change log:
# Version: 1.2.8
# Author: Sander W. van der Laan
# Date: 2026-01-28
# Usage: bash submit_jobs_runPreprocess_rivanna.sh

# Define base directory where studies are stored
# CELLRANGER
# STUDY_DIR="/sfs/gpfs/tardis/project/cphg-millerlab/MetaPlaq/v2/raw_data/post_alignment_outs/cellranger_outs"
# STARSOLO
STUDY_DIR="/sfs/gpfs/tardis/project/cphg-millerlab/MetaPlaq/v2/raw_data/post_alignment_outs/starsolo_outs"

# List of available studies (can be manually updated if needed)
# CELLRANGER
# STUDIES=("Alsaigh_et_al_2022" "Bashore_et_al_2024" "Cheng_et_al" "Chou_et_al_2021" "Eberhardt_et_al_2023" "Fernandez_et_al_2019" "Jaiswal_et_al" "Pan_et_al_2020" "Qian_et_al" "Turner_et_al_2022" "Vacante_et_al" "Wirka_et_al_2019" "Barcia_et_al_2024")
# STARSOLO
STUDIES=("Katyayani_et_al")

# List of platform libraries corresponding to each study (in the same order)
# choices=['10x','smartseq2','celseq2']
# CELLRANGER
# PLATFORM="10x"
# STARSOLO
PLATFORM="smartseq2"

# Loop over studies and submit array jobs
for STUDY in "${STUDIES[@]}"; do
    STUDY_PATH="${STUDY_DIR}/${STUDY}" # Define full path to the study directory, which --input-dir will point to

    if [ -d "$STUDY_PATH" ]; then
        echo "Submitting SLURM job for study: $STUDY" # Define the study name, which --study will point to

        # Get sample count for SLURM array job
        SAMPLE_COUNT=$(find "$STUDY_PATH" -mindepth 1 -maxdepth 1 -type d | wc -l)

        if [ "$SAMPLE_COUNT" -eq 0 ]; then
            echo "Warning: No samples found for study $STUDY. Skipping..."
            continue
        fi

        # Submit SLURM job with correct array size
        sbatch --array=0-$(($SAMPLE_COUNT - 1)) runPreprocess_rivanna.sh --input-dir "$STUDY_PATH" --study "$STUDY" --platform "$PLATFORM" 
    else
        echo "Warning: Study directory $STUDY_PATH does not exist. Skipping..."
    fi
done
