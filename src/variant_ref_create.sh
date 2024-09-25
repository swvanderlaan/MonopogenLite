#!/bin/bash

# Change log:
# * v1.0.5 2024-09-25: Fixed an issue where the bcftools arguments were not properly given.
# * v1.0.4 2024-09-25: Fixed an issue where the AF field was not calculated and filtering was not applied. Added option to choose the AF field. Fixed an issue where the AF was dynamically printed in the file name.
# * v1.0.3 2024-09-25: Fixed an issue where the multi-allelic variants aren't filtered out.
# * v1.0.2 2024-09-24: Added script to remove homozygous calls for chromosome X.
# * v1.0.1 2024-09-24: Added a filter for variants that are homozygous in the phased high-coverage 1000 Genomes VCF files, only for chromosome X.
# * v1.0.0 2024-09-19: Initial version. 
# Version and license information 
VERSION_NAME='Variant Reference Creator'
VERSION='1.0.5'
VERSION_DATE='2024-09-25'
COPYRIGHT='Copyright 1979-2024. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science'
COPYRIGHT_TEXT='''
The MIT License (MIT).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
'''

# Argument parsing function
print_help() {
    echo "$VERSION_NAME version $VERSION ($VERSION_DATE)"
    echo ""
    echo "Usage: $0 --resource-dir <DIR> [--af-field <AF, AF_EUR, AF_EAS, AF_SAS, AF_AMR, AF_AFR>] [--af <FLOAT>] [--variant-type <TYPE>] [--verbose] [--help] [--version]"
    echo ""
    echo "Description:"
    echo "  This script downloads the phased high-coverage 1000 Genomes VCF files "
    echo "  for chromosomes 1-22 and X for use as a variant reference panel."
    echo "  It filters, concatenates, and annotates the VCF files using bcftools."
    echo "  Only bi-allelic SNPs are retained. The output will be stored in the resource "
    echo "  directory."
    echo "  It will run jobs on a SLURM-based system, processing chromosomes 1-22 and X."
    echo "  The script requires the monopogen-environment, and refgenie and bcftools software."
    echo "" 
    echo "  Note that homozygous calls - e.g. 0 or 1 - are not filtered out for chromosome X. "
    echo "  Use the submit_makediploidmalesX.sh to retain only variants with heterozygous calls, "
    echo "  i.e. 0|0, 0|1, and 1|1 are retained."
    echo ""
    echo "Arguments:"
    echo "  --resource-dir  Directory to store resources and outputs."
    echo "  --af-field      Allele frequency field (default: AF)."
    echo "  --af            Allele frequency filter, c.q. variants to keep above AF>X (default: 0.0005)."
    echo "  --variant-type  Variant type filter (default: snp)."
    echo "  --verbose       Enable verbose output."
    echo "  --help          Show this help message and exit."
    echo "  --version       Display version information and exit."
    echo ""
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
    exit 1
}

print_version() {
    echo "$VERSION_NAME version $VERSION ($VERSION_DATE)"
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
    exit 1
}

# Default values
VARIANT_TYPE="snp"
# Assuming AF is an argument for allele frequency field
AF_FIELD="AF"  # Default allele frequency field
AF=0.0005 # 0.5% allele frequency
# Convert AF to scientific notation (e.g., 0.0005 becomes 5e4)
AF_SCI=$(printf "%.0e" $AF)
VERBOSE=0

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --resource-dir) RESOURCE_DIR="$2"; shift ;; # Resource directory
        --af) AF="$2"; AF_SCI=$(printf "%.0e" $AF); shift ;;  # Automatically format AF into scientific notation
        --af-field) AF_FIELD="$2"; shift ;;  # Allow choosing which AF field to use
        --variant-type) VARIANT_TYPE="$2"; shift ;; # Variant type filter
        --verbose) VERBOSE=1 ;;
        --help) print_help ;;
        --version) print_version ;;
        *) echo "Unknown parameter passed: $1"; print_help ;;
    esac
    shift
done

# Starting script
echo "$VERSION_NAME"
echo "version $VERSION ($VERSION_DATE)"
echo ""

# Check if conda is installed
echo "Actvating conda environment..."
source ~/.bashrc
mamba activate monopogen
echo ""

# Check if resource directory is provided
if [[ -z "$RESOURCE_DIR" ]]; then
    echo "Error: --resource-dir flag is required."
    exit 1
fi

# Create the resource directory and output directory
mkdir -p "$RESOURCE_DIR"
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Resource directory created at $RESOURCE_DIR"
fi
OUT_DIR="$RESOURCE_DIR"
mkdir -p "$OUT_DIR"
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Output directory created at $OUT_DIR"
fi

# Set up SLURM parameters
SBATCH_MAIL="FAIL"
SBATCH_MAIL_USER="s.w.vanderlaan-2@umcutrecht.nl"

# Setting up the reference genome
# GRCh38=$(refgenie seek hg38/fasta)
GRCh38="/hpc/dhl_ec/data/references/fasta/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

echo "Starting $VERSION_NAME"
echo ""
echo "These are the settings:"
echo "  Resource directory........: $RESOURCE_DIR"
echo "  Allele frequency field....: $AF_FIELD"
echo "  Allele frequency filter...: $AF ($AF_SCI)"
echo "  Variant type filter.......: $VARIANT_TYPE"
echo "  Reference genome..........: $GRCh38"
echo "  Resource directory........: $RESOURCE_DIR"
echo ""
echo "  Job mail type.............: $SBATCH_MAIL"
echo "  Job mail user.............: $SBATCH_MAIL_USER"
echo ""
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""

# # Submit a job to download data for chromosomes 1-22 and X
# if [[ $VERBOSE -eq 1 ]]; then
#     echo "> Submitting job to download the phased high-coverage 1000 Genomes VCF files."
# fi
# SLURM_DOWNLOAD=$(sbatch --array=1-23 --job-name=pp_download_vcf --output="$RESOURCE_DIR/pp_download_vcf_%A_%a.out" --error="$RESOURCE_DIR/pp_download_vcf_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
# #!/bin/bash
# CHR=\${SLURM_ARRAY_TASK_ID}
# if [[ "\$CHR" == "23" ]]; then
#     wget -P "$RESOURCE_DIR" http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
# else
#     wget -P "$RESOURCE_DIR" http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr\${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
# fi
# EOF
# )

# # Extract the job ID from the output
# SLURM_DOWNLOAD_JOBID=$(echo $SLURM_DOWNLOAD | awk '{print $4}')

# Submit array job to filter chromosomes 1-22 and X
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Submitting job to filter the VCF files."
fi
# SLURM_CREATE=$(sbatch --dependency=afterok:$SLURM_DOWNLOAD_JOBID --array=1-23 --job-name=pp_create --output="$RESOURCE_DIR/pp_create_%A_%a.out" --error="$RESOURCE_DIR/pp_create_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
SLURM_CREATE=$(sbatch --array=1-23 --job-name=pp_create --output="$RESOURCE_DIR/pp_create_%A_%a.out" --error="$RESOURCE_DIR/pp_create_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
source ~/.bashrc
mamba activate monopogen
CHROM=\${SLURM_ARRAY_TASK_ID}

if [[ "\$CHROM" == "23" ]]; then
    echo "Processing chromosome X"
    VCF_IN="${RESOURCE_DIR}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.chrX.vcf.gz"

    echo "> Indexing \$VCF_IN" 
    tabix -fp vcf \$VCF_IN

    echo "> Applying filtering and bi-allelic SNP selection"
    # Step 1: Ensure AF tag is present
    bcftools +fill-tags \$VCF_IN -Oz -o ${OUT_DIR}/temp.chrX.vcf -- -t AF 
    tabix -fp vcf ${OUT_DIR}/temp.chrX.vcf.gz

    # Step 2: Apply filtering and bi-allelic SNP selection
    bcftools view --include '$AF_FIELD>$AF & TYPE="$VARIANT_TYPE"' -m2 -M2 --types snps --regions chrX --output-type z -o ${OUT_DIR}/temp.chrX.vcf \$OUT_FILE
    rm -v ${OUT_DIR}/temp.chrX.vcf.gz
    rm -v ${OUT_DIR}/temp.chrX.vcf.gz.tbi
    
    echo "> Indexing \$OUT_FILE"
    tabix -fp vcf \$OUT_FILE
else
    echo "Processing chromosome \$CHROM"
    VCF_IN="${RESOURCE_DIR}/1kGP_high_coverage_Illumina.chr\${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.chr\${CHROM}.vcf.gz"

    echo "> Indexing \$VCF_IN"
    tabix -fp vcf \$VCF_IN

    echo "> Applying filtering and bi-allelic SNP selection"
    # Step 1: Ensure AF tag is present
    bcftools +fill-tags \$VCF_IN -Oz -o ${OUT_DIR}/temp.chr\${CHROM}.vcf -- -t AF 
    tabix -fp vcf ${OUT_DIR}/temp.chr\${CHROM}.vcf.gz
    
    # Step 2: Apply filtering and bi-allelic SNP selection
    bcftools view --include '$AF_FIELD>$AF & TYPE="$VARIANT_TYPE"' -m2 -M2 --types snps --regions chr\$CHROM --output-type z -o ${OUT_DIR}/temp.\$CHROM.vcf \$OUT_FILE
    rm -v ${OUT_DIR}/temp.chr\${CHROM}.vcf.gz
    rm -v ${OUT_DIR}/temp.chr\${CHROM}.vcf.gz.tbi

    echo "> Indexing \$OUT_FILE"
    tabix -fp vcf \$OUT_FILE
fi
EOF
)

# Extract the job ID from the output
SLURM_CREATE_JOBID=$(echo $SLURM_CREATE | awk '{print $4}')

# Wait for the array job to finish before normalizing the files
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Submitting job to normalize the filtered VCF files."
fi
# SLURM_NORM=$(sbatch --array=1-23 --job-name=pp_norm --output="$RESOURCE_DIR/pp_norm_%A_%a.out" --error="$RESOURCE_DIR/pp_norm_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
SLURM_NORM=$(sbatch --dependency=afterok:$SLURM_CREATE_JOBID --array=1-23 --job-name=pp_norm --output="$RESOURCE_DIR/pp_norm_%A_%a.out" --error="$RESOURCE_DIR/pp_norm_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
source ~/.bashrc
mamba activate monopogen
CHROM=\${SLURM_ARRAY_TASK_ID}
if [[ "\$CHROM" == "23" ]]; then
    CHROM="X"
fi
echo "Normalizing chromosome \$CHROM"
IN_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.chr\${CHROM}.vcf.gz"
OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.chr\${CHROM}.vcf.gz"

echo "> Normalizing \$IN_FILE"
bcftools norm --rm-dup both --check-ref wx --fasta-ref ${GRCh38} \$IN_FILE --output-type z -o \$OUT_FILE

echo "> Indexing \$OUT_FILE"
tabix -fp vcf \$OUT_FILE
EOF
)

# Extract the job ID from the output
SLURM_NORM_JOBID=$(echo $SLURM_NORM | awk '{print $4}')

# Submit the annotation job after concatenation
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Submitting job to annotate the normalized data."
fi
# SLURM_ANNOTATE=$(sbatch --array=1-23 --job-name=pp_annotate --output="$RESOURCE_DIR/pp_annotate_%A_%a.out" --error="$RESOURCE_DIR/pp_annotate_%A_%a.err" --ntasks=1 --cpus-per-task=8 --mem=8G --time=03:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
SLURM_ANNOTATE=$(sbatch --dependency=afterok:$SLURM_NORM_JOBID --array=1-23 --job-name=pp_annotate --output="$RESOURCE_DIR/pp_annotate_%A_%a.out" --error="$RESOURCE_DIR/pp_annotate_%A_%a.err" --ntasks=1 --cpus-per-task=8 --mem=8G --time=03:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
source ~/.bashrc
mamba activate monopogen

# setting chromosome 1-22 and X
CHROM=\${SLURM_ARRAY_TASK_ID}
if [[ "\$CHROM" == "23" ]]; then
    CHROM="X"
fi

echo "Annotating chromosome \$CHROM"
# setting input and output files
IN_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.chr\${CHROM}.vcf.gz"
OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.fixvariantid.chr\${CHROM}.vcf.gz"

# annotate the VCF files
echo "> Annotating \$IN_FILE"
bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' \$IN_FILE --output-type z -o \$OUT_FILE

echo "> Indexing \$OUT_FILE"
tabix -fp vcf \$OUT_FILE
EOF
)

# Extract the job ID from the output
SLURM_ANNOTATE_JOBID=$(echo $SLURM_ANNOTATE | awk '{print $4}')

# Submit the annotation job after concatenation
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Submitting job to subset the normalized and annotated data to list only the variants (for cellsnp)."
fi
# SLURM_SUBSET_CELLSNP=$(sbatch --array=1-23 --job-name=pp_subset --output="$RESOURCE_DIR/pp_subset_%A_%a.out" --error="$RESOURCE_DIR/pp_subset_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
SLURM_SUBSET_CELLSNP=$(sbatch --dependency=afterok:$SLURM_ANNOTATE_JOBID --array=1-23 --job-name=pp_subset --output="$RESOURCE_DIR/pp_subset_%A_%a.out" --error="$RESOURCE_DIR/pp_subset_%A_%a.err" --ntasks=1 --cpus-per-task=8 --mem=8G --time=03:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
source ~/.bashrc
mamba activate monopogen

# setting chromosome 1-22 and X
CHROM=\${SLURM_ARRAY_TASK_ID}
if [[ "\$CHROM" == "23" ]]; then
    CHROM="X"
fi

echo "Subsetting chromosome \$CHROM"
# setting input and output files
IN_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.fixvariantid.chr\${CHROM}.vcf.gz"
OUT_FILE_CELLSNP="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.fixvariantid.chr\${CHROM}.cellsnp.vcf.gz"

# subsetting the VCF files
echo "> Subsetting \$IN_FILE"
bcftools view --samples "." \$IN_FILE --output-type z -o \$OUT_FILE_CELLSNP --force-samples

echo "> Indexing \$OUT_FILE_CELLSNP"
tabix -fp vcf \$OUT_FILE_CELLSNP
EOF
)

# Extract the job ID from the output
SLURM_SUBSET_CELLSNP_JOBID=$(echo $SLURM_SUBSET_CELLSNP | awk '{print $4}')

# Wait for the array job to finish before concatenating the files
if [[ $VERBOSE -eq 1 ]]; then
    echo "> Submitting job to concatenate the filtered VCF files for cellsnp."
fi
# SLURM_CONCAT_CELLSNP=$(sbatch --job-name=pp_concat --output="$RESOURCE_DIR/pp_concat.out" --error="$RESOURCE_DIR/pp_concat.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=03:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
SLURM_CONCAT_CELLSNP=$(sbatch --dependency=afterok:$SLURM_SUBSET_CELLSNP_JOBID --job-name=pp_concat --output="$RESOURCE_DIR/pp_concat.out" --error="$RESOURCE_DIR/pp_concat.err" --ntasks=1 --cpus-per-task=8 --mem=16G --time=02:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
source ~/.bashrc
mamba activate monopogen

echo "Concatenating the filtered VCF files for cellsnp"
# Use shell expansion to capture the files
chrom_files=(${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.fixvariantid.chr*.cellsnp.vcf.gz)

# Check if the files exist
if [[ \${#chrom_files[@]} -eq 0 ]]; then
    echo "No VCF files found matching the pattern."
    exit 1
fi

# Concatenate the files
OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_${AF_FIELD}_${AF_SCI}.norm.fixvariantid.chr1_22X.cellsnp.vcf.gz"

echo "> Concatenating \${chrom_files[@]}"
bcftools concat "\${chrom_files[@]}" --output-type z -o \$OUT_FILE

# Index the output
echo "> Indexing \$OUT_FILE"
tabix -fp vcf \$OUT_FILE
EOF
)

if [[ $VERBOSE -eq 1 ]]; then
    echo ""
    echo "All jobs submitted. Let's wait and see. Have a beer, buddy!"
fi


echo ""
print_version
### END OF SCRIPT ###

