#!/bin/bash

# Change log:
# * v1.0.1 2024-09-30: Fixed references. Fixed array-submission. Added --dry-run, --debug modes.
# * v1.0.0 2024-09-19: Initial version. 
# Version and license information 
VERSION_NAME='Variant Reference Checker'
VERSION='1.0.1'
VERSION_DATE='2024-09-30'
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
    echo "Usage: $0 --input <full-path-to-input.vcf.gz> [--dry-run] [--verbose] [--help] [--version] [--job-name <variant_ref_checker>] [--cpus <1>] [--mem <memory>] [--time <time>] [--mailtype <mail type>] [--mailuser <email>]"
    echo
    echo "  --input         Full path to the input VCF file (including .vcf.gz). Required."
    echo "  --job-name      SLURM job name (default: variant_ref_checker). Optional."
    echo "  --cpus          Number of CPUs per task (default: 1). Optional."
    echo "  --mem           Memory per task (default: 8G). Optional."
    echo "  --time          Time limit per task (default: 01:00:00). Optional."
    echo "  --mailtype      Mail type for SLURM (default: FAIL). Optional."
    echo "  --mailuser      Email for SLURM notifications (default: s.w.vanderlaan-2@umcutrecht.nl). Optional."
    echo "  --dry-run       Perform a dry run without submitting the job. Optional."
    echo "  --verbose       Enable verbose output. Optional."
    echo "  --debug         Enable debug output. Optional."
    echo "  --help          Show this help message and exit."
    echo "  --version       Show the script version and exit."
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
SBATCH_JOB_NAME="variant_ref_checker"
SBATCH_CPUS=1
SBATCH_MEM="8G"
SBATCH_TIME="01:00:00"
SBATCH_MAILTYPE="FAIL"
SBATCH_MAILUSER="s.w.vanderlaan-2@umcutrecht.nl"
VERBOSE=0
DEBUG=0
DRY_RUN=0

# MonopogenLite location
MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"

# Parse command line arguments
echo "> Parsing command line arguments."
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_VCF="$2"; shift ;;
        --job-name) SBATCH_JOB_NAME="$2"; shift ;;
        --cpus) SBATCH_CPUS="$2"; shift ;;
        --mem) SBATCH_MEM="$2"; shift ;;
        --time) SBATCH_TIME="$2"; shift ;;
        --mailtype) SBATCH_MAILTYPE="$2"; shift ;;
        --mailuser) SBATCH_MAILUSER="$2"; shift ;;
        --dry-run) DRY_RUN=1 ;;
        --verbose) VERBOSE=1 ;;
        --debug) DEBUG=1 ;;
        --version) print_version ;;
        --help) print_help ;;
        *) echo "Unknown parameter passed: $1"; print_help ;;
    esac
    shift
done

# Check mandatory parameters
if [[ "$DEBUG" -eq 1 ]]; then
    echo "> Checking mandatory parameters."
fi
if [[ -z "$INPUT_VCF" ]]; then
    echo "Error: --input is required."
    print_help
    exit 1
fi

# Check if input file exists
if [[ "$DEBUG" -eq 1 ]]; then
    echo "> Checking if input file exists."
fi
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "Error: Input file '$INPUT_VCF' does not exist."
    exit 1
fi

# Create the SLURM batch job script
SBATCH_SCRIPT="$MPG/submit_variant_ref_checker.sbatch"

# Extract directory and base file name
if [[ "$DEBUG" -eq 1 ]]; then
    echo "> Extracting the directory and base file name from the input VCF file."
fi
VCF_DIR=$(dirname "$INPUT_VCF")
BASE_NAME=$(basename "$INPUT_VCF")

# Get the part of the file name before ".chr#."
if [[ "$DEBUG" -eq 1 ]]; then
    echo "> Extracting the base prefix from the input VCF file name."
fi
PREFIX=$(echo "$BASE_NAME" | sed -r 's/\.chr[0-9XY]+\.vcf\.gz//')

if [[ -z "$PREFIX" ]]; then
    echo "Error: Unable to parse the input VCF file name to determine the base prefix."
    exit 1
fi

echo "Starting $VERSION_NAME"
echo ""
echo "These are the settings:"
echo "  Input file................: $INPUT_VCF"
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
echo "  Dry run mode..............: $DRY_RUN"
echo "  Debug mode................: $DEBUG"
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""

# Create the SLURM job array
echo "> Creating the SLURM job array script."

cat << EOF > $SBATCH_SCRIPT
#!/bin/bash

#SBATCH --job-name=$SBATCH_JOB_NAME
#SBATCH --array=1-23
#SBATCH --cpus-per-task=$SBATCH_CPUS
#SBATCH --mem=$SBATCH_MEM
#SBATCH --time=$SBATCH_TIME
#SBATCH --mail-type=$SBATCH_MAILTYPE
#SBATCH --mail-user=$SBATCH_MAILUSER
#SBATCH --output=${SBATCH_JOB_NAME}_%A_%a.out
#SBATCH --error=${SBATCH_JOB_NAME}_%A_%a.err

source ~/.bashrc
mamba activate monopogen

echo "$VERSION_NAME"
echo "version $VERSION ($VERSION_DATE)"
echo ""
echo "These are the settings:"
echo "  Input file................: $INPUT_VCF"
echo "  VCF directory.............: $VCF_DIR"
echo "  VCF file prefix...........: $PREFIX"
echo ""
echo "  SLURM job name............: $SBATCH_JOB_NAME"
echo "  SLURM Array ID............: \$SLURM_ARRAY_TASK_ID"
echo "  SLURM CPUs................: $SBATCH_CPUS"
echo "  SLURM memory..............: $SBATCH_MEM"
echo "  SLURM time................: $SBATCH_TIME"
echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
echo "  SLURM mail user...........: $SBATCH_MAILUSER"
echo ""
echo "  Dry run mode..............: $DRY_RUN"
echo "  Debug mode................: $DEBUG"
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""
echo "Running $VERSION_NAME..."

echo "> Chromosome mapping for SLURM array task ID # \$SLURM_ARRAY_TASK_ID..."
CHR_LIST=($(seq 1 22) "X")

echo "> Set the chromosome based on the array task ID..."
CHR=\${CHR_LIST[\$SLURM_ARRAY_TASK_ID - 1]}

echo "  - Construct the input VCF file for the current chromosome..."
CHR_VCF="${VCF_DIR}/${PREFIX}.chr\${CHR}.vcf.gz"

echo "> Check if the file exists..."
if [[ ! -f "\$CHR_VCF" ]]; then
    echo "Error: VCF file not found for chromosome \${CHR}: \$CHR_VCF"
    exit 1
fi

echo "> Run the Python script for the specified chromosome..."
python3 $MPG/src/variant_ref_checker.py --input "\$CHR_VCF" --chr \$CHR $( [[ $VERBOSE -eq 1 ]] && echo "--verbose" )

if [ \$? -eq 0 ]; then
  echo "$VERSION_NAME finished successfully. Let's have a beer, buddy!"
else
  echo "$VERSION_NAME encountered an error."
fi

mamba deactivate

EOF

# Make the script executable
echo "> Make the script executable."
chmod +x $SBATCH_SCRIPT

# Submit the job to SLURM
echo ""
echo "Submitting the job to SLURM..."
if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN: sbatch $SBATCH_SCRIPT"
    exit 0
else 
    JOB_ID=$(sbatch $SBATCH_SCRIPT | awk '{print $4}')
    echo ">> Job submitted with ID: $JOB_ID"
fi

echo ""
if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN completed. No jobs were submitted."
else
    echo "All jobs submitted, this will take a while. Let's have a beer, buddy!"
fi

echo ""
print_version
### END OF SCRIPT ###
