#!/bin/bash

VERSION_NAME='Variant Reference Creator'
VERSION='1.0.0'
VERSION_DATE='2024-09-18'
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
    echo "Usage: $0 --resource-dir <DIR> [--verbose] [--help] [--version]"
    echo ""
    echo "Description:"
    echo "  This script downloads the phased high-coverage 1000 Genomes VCF files "
    echo "  for chromosomes 1-22 and X for use as a variant reference panel."
    echo "  This script filters, concatenates, and annotates the VCF files using bcftools."
    echo "  It runs jobs on a SLURM-based system, processing chromosomes 1-22 and X."
    echo ""
    echo "Arguments:"
    echo "  --resource-dir  Directory to store resources and outputs."
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
    exit 1
}

# Default values
VERBOSE=0

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --resource-dir) RESOURCE_DIR="$2"; shift ;;
        --verbose) VERBOSE=1 ;;
        --help) print_help ;;
        --version) print_version ;;
        *) echo "Unknown parameter passed: $1"; print_help ;;
    esac
    shift
done

# Check if resource directory is provided
if [[ -z "$RESOURCE_DIR" ]]; then
    echo "Error: --resource-dir flag is required."
    exit 1
fi

# Create the resource directory
mkdir -p "$RESOURCE_DIR"
if [[ $VERBOSE -eq 1 ]]; then
    echo "Resource directory created at $RESOURCE_DIR"
fi

SBATCH_MAIL="FAIL"
SBATCH_MAIL_USER="s.w.vanderlaan-2@umcutrecht.nl"

# SLURM parameters
OUT_DIR="$RESOURCE_DIR/output"
mkdir -p "$OUT_DIR"

# Submit a job to download data for chromosomes 1-22 and X
SLURM_DOWNLOAD=$(sbatch --array=1-23 --job-name=pp_download_vcf --output="$RESOURCE_DIR/pp_download_vcf_%A_%a.out" --error="$RESOURCE_DIR/pp_download_vcf_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
CHR=\${SLURM_ARRAY_TASK_ID}
if [[ "\$CHR" == "23" ]]; then
    wget -P "$RESOURCE_DIR" http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
else
    wget -P "$RESOURCE_DIR" http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr\${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
fi
EOF
)

# Extract the job ID from the output
SLURM_DOWNLOAD_JOBID=$(echo $SLURM_DOWNLOAD | awk '{print $4}')

# Submit array job to filter chromosomes 1-22 and X
SLURM_CREATE=$(sbatch --dependency=afterok:$SLURM_DOWNLOAD_JOBID --array=1-23 --job-name=pp_create --output="$RESOURCE_DIR/pp_create_%A_%a.out" --error="$RESOURCE_DIR/pp_create_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
CHROM=\${SLURM_ARRAY_TASK_ID}
if [[ "\$CHROM" == "23" ]]; then
    VCF_IN="${RESOURCE_DIR}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af.chrX.vcf.gz"
    bcftools view -i 'AF>0.0 & TYPE="snp"' -s "." \$VCF_IN -r chrX -Oz -o \$OUT_FILE --force-samples
else
    VCF_IN="${RESOURCE_DIR}/1kGP_high_coverage_Illumina.chr\${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af.chr\${CHROM}.vcf.gz"
    bcftools view -i 'AF>0.0 & TYPE="snp"' -s "." \$VCF_IN -r chr\$CHROM -Oz -o \$OUT_FILE --force-samples
fi
EOF
)

# Extract the job ID from the output
SLURM_CREATE_JOBID=$(echo $SLURM_CREATE | awk '{print $4}')

# Wait for the array job to finish before concatenating the files
SLURM_CONCAT=$(sbatch --dependency=afterok:$SLURM_CREATE_JOBID --job-name=pp_concat --output="$RESOURCE_DIR/pp_concat.out" --error="$RESOURCE_DIR/pp_concat.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af.chr1_22X.vcf.gz"
chrom_files=\$(ls "${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af.chr*.vcf.gz" | tr '\n' ' ')
bcftools concat \$chrom_files -Oz -o \$OUT_FILE
EOF
)

# Extract the job ID from the output
SLURM_CONCAT_JOBID=$(echo $SLURM_CONCAT | awk '{print $4}')

# Submit the annotation job after concatenation
SLURM_ANNOTATE=$(sbatch --dependency=afterok:$SLURM_CONCAT_JOBID --job-name=pp_annotate --output="$RESOURCE_DIR/pp_annotate.out" --error="$RESOURCE_DIR/pp_annotate.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
INPUT="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af.chr1_22X.vcf.gz"
OUT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af_5e4.chr1_22X.vcf.gz"
bcftools annotate -x ^INFO/AF -i 'AF>0.0005' \$INPUT -Oz -o \$OUT_FILE
EOF
)

# Extract the job ID from the output
SLURM_ANNOTATE_JOBID=$(echo $SLURM_ANNOTATE | awk '{print $4}')

# Submit job to split the concatenated VCF back into individual chromosomes (1-22 and X as 23)
sbatch --dependency=afterok:$SLURM_ANNOTATE_JOBID --array=1-23 --job-name=pp_split --output="$RESOURCE_DIR/pp_split_%A_%a.out" --error="$RESOURCE_DIR/pp_split_%A_%a.err" --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00 --mail-type="$SBATCH_MAIL" --mail-user="$SBATCH_MAIL_USER" << EOF
#!/bin/bash
CHR=\${SLURM_ARRAY_TASK_ID}
if [[ "\$CHR" == "23" ]]; then
    CHR="X"
fi
OUT_SPLIT_FILE="${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af_5e4.chr\${CHR}.vcf.gz"
bcftools view -r \$CHR ${OUT_DIR}/1kGP_high_coverage_Illumina.SNVonly_poly.filtered_af_5e4.chr1_22X.vcf.gz -Oz -o \$OUT_SPLIT_FILE
EOF

if [[ $VERBOSE -eq 1 ]]; then
    echo "Jobs submitted."
fi
