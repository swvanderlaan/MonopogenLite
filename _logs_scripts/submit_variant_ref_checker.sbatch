#!/bin/bash

#SBATCH --job-name=variant_ref_checker
#SBATCH --array=1-23
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl
#SBATCH --output=variant_ref_checker_%A_%a.out
#SBATCH --error=variant_ref_checker_%A_%a.err

source ~/.bashrc
mamba activate monopogen

echo "Variant Reference Checker"
echo "version 1.0.2 (2024-10-01)"
echo ""
echo "These are the settings:"
echo "  Input file................: /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.chr1.vcf.gz"
echo "  VCF directory.............: /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources"
echo "  VCF file prefix...........: 1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid"
echo ""
echo "  SLURM job name............: variant_ref_checker"
echo "  SLURM Array ID............: $SLURM_ARRAY_TASK_ID"
echo "  SLURM CPUs................: 1"
echo "  SLURM memory..............: 8G"
echo "  SLURM time................: 01:00:00"
echo "  SLURM mail type...........: FAIL"
echo "  SLURM mail user...........: s.w.vanderlaan-2@umcutrecht.nl"
echo ""
echo "  Dry run mode..............: 0"
echo "  Debug mode................: 0"
echo "  Verbosity.................: 1"
echo "  Version...................: 1.0.2 (2024-10-01)"
echo ""
echo "Running Variant Reference Checker..."

echo "> Chromosome mapping for SLURM array task ID # [$SLURM_ARRAY_TASK_ID]..."
CHR_LIST=(1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22 "X")

echo "> Set the chromosome based on the array task ID..."
CHR=${CHR_LIST[$SLURM_ARRAY_TASK_ID - 1]}

echo "  - Construct the input VCF file for the current chromosome # [$CHR]..."
CHR_VCF="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.chr${CHR}.vcf.gz"

echo "> Check if the file exists..."
if [[ ! -f "$CHR_VCF" ]]; then
    echo "Error: VCF file not found for chromosome ${CHR}: $CHR_VCF"
    exit 1
fi

echo "> Run the Python script for the specified chromosome..."
python3 /hpc/local/Rocky8/dhl_ec/software/MonopogenLite/src/variant_ref_checker.py --input "$CHR_VCF" --chr $CHR --verbose

if [ $? -eq 0 ]; then
  echo "Variant Reference Checker finished successfully. Let's have a beer, buddy!"
else
  echo "Variant Reference Checker encountered an error."
fi

mamba deactivate

