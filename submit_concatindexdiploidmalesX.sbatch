#!/bin/bash
#SBATCH --job-name=concatdiploidmalesX
#SBATCH --partition=standard
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl
#SBATCH --output=concatdiploidmalesX_%j.out
#SBATCH --error=concatdiploidmalesX_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=02:00:00

source ~/.bashrc
source ~/.bash_profile
micromamba activate monopogen

DEBUG_FLAG=0

echo "Submit MakeDiploidMalesX"
echo "version 1.1.8 (2024-10-30)"
echo ""
echo "These are the settings:"
echo "  Input file................: /project/cphg-millerlab/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.chrX.vcf.gz"
echo "  Output file...............: /project/cphg-millerlab/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.makediploid.chrX.vcf.gz"
echo "  Chunk size................: 9000"
if [[ -n "" ]]; then
    echo "  Changes file..............: "
fi
if [[ -n "" ]]; then
    echo "  Reverse file..............: "
fi
echo ""
echo "  SLURM CPUs................: 1"
echo "  SLURM partition...........: standard"
echo "  SLURM time................: 12:00:00"
if [[ standard == "gpu" ]]; then
    echo "  SLURM mem GPU.............: 16G"
    echo "  SLURM GPUs per node.......: 1"
    echo "  SLURM tasks...............: 1"
    echo "  SLURM account.............: cphg-millerlab"
else 
    echo "  SLURM memory..............: 16G"
fi
echo "  SLURM mail type...........: FAIL"
echo "  SLURM mail user...........: s.w.vanderlaan-2@umcutrecht.nl"
echo ""
echo "  Dry run mode..............: 0"
echo "  Debug mode................: 0"
echo "  Verbosity.................: 0"
echo "  Version...................: 1.1.8 (2024-10-30)"
echo ""
echo "Concatenating processed VCF files."

echo "> Concatenate the processed VCF files..."
echo "  - Construct the ordered list of files."
VCF_LIST=""
CHUNKSIZE=9000
for CHUNK in $(seq 1 $CHUNKSIZE); do
    FILE="/project/cphg-millerlab/software/MonopogenLite/resources/chrX.part${CHUNK}_processed.vcf.gz"
    echo "  - Processing chunk $CHUNK and adding [$FILE] to the list..."
    if [[ -f "$FILE" ]]; then
        VCF_LIST="$VCF_LIST $FILE"
    else
        echo "ERROR: File $FILE does not exist..."
        exit 1
    fi
done

if [[ $DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Listed all the processed VCF files:"
    printf "%s\n" "$VCF_LIST"
fi

bcftools concat -Oz -o /project/cphg-millerlab/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.makediploid.chrX.vcf.gz $VCF_LIST
if [ 0 -ne 0 ]; then
    echo "ERROR: Concatenation failed."
    exit 1
fi

echo "> Index the concatenated VCF file..."
tabix -fp vcf /project/cphg-millerlab/software/MonopogenLite/resources/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.makediploid.chrX.vcf.gz
if [ 0 -ne 0 ]; then
    echo "ERROR: Indexing with tabix failed."
    exit 1
fi

echo "> Removing the processed VCF files..."
for CHUNK in $(seq 1 $CHUNKSIZE); do
    FILE="/project/cphg-millerlab/software/MonopogenLite/resources/chrX.part${CHUNK}_processed.vcf.gz"
    FILE_TBI="/project/cphg-millerlab/software/MonopogenLite/resources/chrX.part${CHUNK}_processed.vcf.gz.tbi"
    FILE_RAW="/project/cphg-millerlab/software/MonopogenLite/resources/chrX.part${CHUNK}.vcf.gz"
    FILE_RAW_TBI="/project/cphg-millerlab/software/MonopogenLite/resources/chrX.part${CHUNK}.vcf.gz.tbi"
    echo "  - Removing processed VCF file: [$FILE]"
    rm -v $FILE
    echo "  - Removing processed VCF file: [$FILE_TBI]"
    rm -v $FILE_TBI
    echo "  - Removing raw VCF file: [$FILE_RAW]"
    rm -v $FILE_RAW
    echo "  - Removing raw VCF file: [$FILE_RAW_TBI]"
    rm -v $FILE_RAW_TBI
done

if [ $? -eq 0 ]; then
    echo ""
    echo "Submit MakeDiploidMalesX finished successfully. Let's have a beer, buddy!"
    echo ""
else
    echo ""
    echo "ERROR: Submit MakeDiploidMalesX encountered an error."
    echo ""
fi

micromamba deactivate

