#!/bin/bash

# Change log:
# * v1.0.0 2024-09-19: Initial version. 
# Version and license information 
VERSION_NAME='Variant Reference Creator'
VERSION='1.0.0'
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
    echo "Usage: $0 --input <full-path-to-vcf-file> [--verbose] [--help] [--version] [--mem <memory>] [--time <time>] [--mailtype <mail type>] [--mailuser <email>]"
    echo
    echo "  --input         Full path to the input VCF file (including .vcf.gz)."
    echo "  --verbose       Enable verbose output (optional)."
    echo "  --help          Show this help message and exit."
    echo "  --version       Show the script version and exit."
    echo "  --mem           Memory per task (default: 4G)."
    echo "  --time          Time limit per task (default: 00:30:00)."
    echo "  --mailtype      Mail type for SLURM (default: FAIL)."
    echo "  --mailuser      Email for SLURM notifications (default: s.w.vanderlaan-2@umcutrecht.nl)."
    exit 0
}

print_version() {
    echo "$VERSION_NAME version $VERSION ($VERSION_DATE)"
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
    exit 0
}

# Starting script
echo "$VERSION_NAME"
echo "version $VERSION ($VERSION_DATE)"
echo ""

# Check if conda is installed
echo "Activating conda environment..."
source ~/.bashrc
mamba activate monopogen
echo ""

# Default values for optional parameters
SBATCH_JOB_NAME="makediploidmalesX"
SBATCH_CPUS=1
SBATCH_MEM="4G"
SBATCH_TIME="00:30:00"
SBATCH_MAILTYPE="FAIL"
SBATCH_MAILUSER="s.w.vanderlaan-2@umcutrecht.nl"
VERBOSE=0

# MonopogenLite location
MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_VCF="$2"; shift ;;
        --job-name) SBATCH_JOB_NAME="$2"; shift ;;
        --cpus) SBATCH_CPUS="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --time) TIME="$2"; shift ;;
        --mailtype) MAILTYPE="$2"; shift ;;
        --mailuser) MAILUSER="$2"; shift ;;
        --verbose) VERBOSE=1 ;;
        --version) print_version ;;
        --help) print_help ;;
        *) echo "Unknown parameter passed: $1"; print_help ;;
    esac
    shift
done

# Check mandatory parameters
if [[ -z "$INPUT_VCF" ]]; then
    echo "Error: --input is required."
    usage
    exit 1
fi

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

# Create the SLURM batch job script
SBATCH_SCRIPT="$MPG/submit_makediploidmalesX.sbatch"

echo "Starting $VERSION_NAME"
echo ""
echo "These are the settings:"
echo "  Input file................: $INPUT_FILE"
echo "  SLURM job name............: $SBATCH_JOB_NAME"
echo "  SLURM CPUs................: $SBATCH_CPUS"
echo "  SLURM memory..............: $SBATCH_MEM"
echo "  SLURM time................: $SBATCH_TIME"
echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
echo "  SLURM mail user...........: $SBATCH_MAILUSER"
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""

# Extract directory and base file name
VCF_DIR=$(dirname "$INPUT_VCF")
BASE_NAME=$(basename "$INPUT_VCF")

# Get the part of the file name before ".chr#."
PREFIX=$(echo "$BASE_NAME" | sed -r 's/\.chr[0-9XY]+\.vcf\.gz//')

if [[ -z "$PREFIX" ]]; then
    echo "Error: Unable to parse the input VCF file name to determine the base prefix."
    exit 1
fi

# Submit the SLURM job array
cat << EOF > $SBATCH_SCRIPT
#!/bin/bash

#SBATCH --job-name=$SBATCH_JOB_NAME
#SBATCH --cpus-per-task=$SBATCH_CPUS
#SBATCH --mem=$SBATCH_MEM
#SBATCH --time=$SBATCH_TIME
#SBATCH --mail-type=$SBATCH_MAILTYPE
#SBATCH --mail-user=$SBATCH_MAILUSER
#SBATCH --output=${SBATCH_JOB_NAME}_%A_%a.out
#SBATCH --error=${SBATCH_JOB_NAME}_%A_%a.err
#SBATCH --array=1-23

source ~/.bashrc
mamba activate monopogen

echo "$VERSION_NAME"
echo "version $VERSION ($VERSION_DATE)"
echo ""
echo "These are the settings:"
echo "  Input file................: $INPUT_FILE"
echo "  VCF directory.............: $VCF_DIR"
echo "  VCF file prefix...........: $PREFIX"
echo ""
echo "  SLURM job name............: $SBATCH_JOB_NAME"
echo "  SLURM CPUs................: $SBATCH_CPUS"
echo "  SLURM memory..............: $SBATCH_MEM"
echo "  SLURM time................: $SBATCH_TIME"
echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
echo "  SLURM mail user...........: $SBATCH_MAILUSER"
echo ""
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""
echo "Running $VERSION_NAME..."

echo "> Chromosome mapping for SLURM array task ID 23..."
CHR_LIST=($(seq 1 22) "X")

echo "> Set the chromosome based on the array task ID..."
CHR=\${CHR_LIST[\$SLURM_ARRAY_TASK_ID - 1]}

echo "Construct the input VCF file for the current chromosome..."
CHR_VCF="${VCF_DIR}/${PREFIX}.chr\${CHR}.vcf.gz"

echo "> Check if the file exists..."
if [[ ! -f "\$CHR_VCF" ]]; then
    echo "Error: VCF file not found for chromosome \${CHR}: \$CHR_VCF"
    exit 1
fi

echo "> Run the Python script for the specified chromosome..."
python3 bcftools_stats_plot.py --input "\$CHR_VCF" ${VERBOSE}

if [ \$? -eq 0 ]; then
  echo "$VERSION_NAME finished successfully. Let's have a beer, buddy!"
else
  echo "$VERSION_NAME encountered an error."
fi

mamba deactivate

EOF

# Make the script executable
chmod +x $SBATCH_SCRIPT

# Submit the job to SLURM
JOB_ID=$(sbatch $SBATCH_SCRIPT | awk '{print $4}')
echo "Job submitted with ID: $JOB_ID"

echo ""
print_version
### END OF SCRIPT ###