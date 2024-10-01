#!/bin/bash
#SBATCH --job-name=makediploidmalesX
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl
#SBATCH --output=makediploidmalesX_%A_%a.out
#SBATCH --error=makediploidmalesX_%A_%a.err

source ~/.bashrc
mamba activate monopogen

DEBUG_FLAG=1

echo "Submit MakeDiploidMalesX"
echo "version 1.1.6 (2024-09-30)"
echo ""
echo "These are the settings:"
echo "  Input file................: /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.chrX.vcf.gz"
echo "  Output file...............: /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.makediploid.chrX.vcf.gz"
echo "  Chunk size................: 100"
if [[ -n "" ]]; then
    echo "  Changes file..............: "
fi
if [[ -n "" ]]; then
    echo "  Reverse file..............: "
fi
echo ""
echo "  SLURM CPUs................: 4"
echo "  SLURM Array ID............: $SLURM_ARRAY_TASK_ID (of 100)"
echo "  SLURM memory..............: 16G"
echo "  SLURM time................: 01:00:00"
echo "  SLURM mail type...........: FAIL"
echo "  SLURM mail user...........: s.w.vanderlaan-2@umcutrecht.nl"
echo ""
echo "  Dry run mode..............: 0"
echo "  Debug mode................: 1"
echo "  Verbosity.................: 0"
echo "  Version...................: 1.1.6 (2024-09-30)"
echo ""
echo "Making converting haploid genotypes to diploid genotypes."

echo "> Extract the region for this SLURM task (# $SLURM_ARRAY_TASK_ID) from the BED file."

echo "  - Get the chunk number from the SLURM array task ID and the bed file..."
CHUNK_NUMBER_RAW=$SLURM_ARRAY_TASK_ID # SLURM array task ID
CHUNK_NUMBER=${CHUNK_NUMBER_RAW}p  # adding p to the number
BED_FILE="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX_100pieces.bed" # BED file with chunks

# Debugs
if [[ $DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Extracting chunk number: $CHUNK_NUMBER from $CHUNK_NUMBER_RAW"
fi
if [[ $DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Extracting region from BED file: $BED_FILE"
fi

echo "  - Extract the region from the BED file..."
RAW_REGION=$(sed -n "$CHUNK_NUMBER" "$BED_FILE")
RAW_REGION=$(echo "$RAW_REGION" | sed 's/^ *//;s/ *$//')  # Clean up leading/trailing whitespace
if [[ $DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Extracted raw region: '$RAW_REGION' (chunk #: $CHUNK_NUMBER of 100 for job $CHUNK_NUMBER_RAW)"
fi
echo "  - Check if a valid region was extracted -- should be '$RAW_REGION'..."
if [[ -z "$RAW_REGION" ]]; then
  echo "ERROR: No valid region extracted for chunk # $SLURM_ARRAY_TASK_ID (or $CHUNK_NUMBER_RAW) from '$BED_FILE'."
  exit 1
fi

echo "  - Prepare the region for bcftools..."
REGION=$(echo "$RAW_REGION" | awk '{print $1 ":" $2 "-" $3}')
echo "  - Check if the region ($REGION) was successfully extracted..."
if [[ -z "$REGION" ]]; then
  echo "ERROR: RAW_REGION was: '$RAW_REGION'"
  echo "ERROR: No region found for SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID."
  exit 1
fi

echo "> Processing region '$REGION'..."

echo "  - Process the specific region (chunk# $CHUNK_NUMBER_RAW) using bcftools view..."
bcftools view -r $REGION /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.chrX.vcf.gz -Oz -o /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX.part${CHUNK_NUMBER_RAW}.vcf.gz

echo "  - Index the VCF file..."
tabix -fp vcf /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX.part${CHUNK_NUMBER_RAW}.vcf.gz

# Debugs
if [[ $DEBUG_FLAG -eq "1" ]]; then
    echo "DEBUG: Check if the VCF file was successfully processed..."
    bcftools view /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX.part${CHUNK_NUMBER_RAW}.vcf.gz | head
    bcftools stats /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX.part${CHUNK_NUMBER_RAW}.vcf.gz
fi

echo "  - Fix haploid genotypes in males..."
python3 /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/src/makediploidmalesX.py --input-file /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX.part${CHUNK_NUMBER_RAW}.vcf.gz --output-file /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/chrX.part${CHUNK_NUMBER_RAW}_processed.vcf.gz   

if [ $? -eq 0 ]; then
    echo ""
    echo "Submit MakeDiploidMalesX finished successfully. Let's have a beer, buddy!"
    echo ""
else
    echo ""
    echo "ERROR: Submit MakeDiploidMalesX encountered an error."
    echo ""
fi

mamba deactivate

