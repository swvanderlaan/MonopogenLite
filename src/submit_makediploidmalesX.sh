#!/bin/bash

# Change log:
# * v1.0.5 2024-09-27: Changed default values for SLURM.
# * v1.0.4 2024-09-25: Added a check if the input file exists. Added optional --changes and --reverse flags.
# * v1.0.3 2024-09-25: Changed script name.
# * v1.0.2 2024-09-24: Fixed issue where there was no --version flag in the help message.
# * v1.0.1 2024-09-24: Added a filter for variants that are homozygous in the phased high-coverage 1000 Genomes VCF files, only for chromosome X.
# * v1.0.0 2024-09-24: Initial version. 
# Version and license information 
VERSION_NAME='Submit MakeDiploidMalesX'
VERSION='1.0.5'
VERSION_DATE='2024-09-27'
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
    echo "Usage: $0 --input <input.vcf.gz> --output <output.vcf.gz> [--changes <changes.txt.gz>] [--reverse <reverse.txt.gz>] [--job-name <job_name>] [--cpus <num_cpus>] [--mem <memory>] [--time <time>] [--mail <mail-type>] [--user <mail-user>] [--verbose]"
    echo ""
    echo "Description:"
    echo "  This script will submit a job to the SLURM scheduler to make haploid genotypes in males diploid given a VCF file."
    echo ""
    echo "Arguments:"
    echo "  --input         The input VCF file."
    echo "  --output        The output VCF file."
    echo "  --changes       Gzipped file to save the list of changes (optional)."
    echo "  --reverse       Gzipped file with list of changes to reverse (optional)."
    echo "  --job-name      The name of the SLURM job."
    echo "  --cpus          The number of CPUs to use."
    echo "  --mem           The amount of memory to use."
    echo "  --time          The maximum time to run the job."
    echo "  --mail          The type of mail to send."
    echo "  --user          The email address to send the mail to."
    echo "  --verbose       Enable verbose output."
    echo "  --help          Show this help message and exit."
    echo "  --version       Display version information and exit."
    echo ""
    echo "$COPYRIGHT"
    echo "$COPYRIGHT_TEXT"
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

# Default values
SBATCH_JOB_NAME="makediploidmalesX"
SBATCH_CPUS=4
SBATCH_MEM="32G"
SBATCH_TIME="12:00:00"
SBATCH_MAILTYPE="FAIL"
SBATCH_MAILUSER="s.w.vanderlaan-2@umcutrecht.nl"
VERBOSE=0
CHANGES_FILE=""
REVERSE_FILE=""

# MonopogenLite location
MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_FILE="$2"; shift ;;
        --output) OUTPUT_FILE="$2"; shift ;;
        --changes) CHANGES_FILE="$2"; shift ;;
        --reverse) REVERSE_FILE="$2"; shift ;;
        --job-name) SBATCH_JOB_NAME="$2"; shift ;;
        --cpus) SBATCH_CPUS="$2"; shift ;;
        --mem) SBATCH_MEM="$2"; shift ;;
        --time) SBATCH_TIME="$2"; shift ;;
        --mail) SBATCH_MAILTYPE="$2"; shift ;;
        --user) SBATCH_MAILUSER="$2"; shift ;;
        --verbose) VERBOSE=1 ;;
        --version) print_version ;;
        --help) print_help ;;
        *) echo "Unknown parameter passed: $1"; print_help ;;
    esac
    shift
done

# Check if input and output are provided
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
    echo "Error: Both --input and --output arguments are required."
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
echo "  Output file...............: $OUTPUT_FILE"
echo "  Changes file..............: $CHANGES_FILE"
echo "  Reverse file..............: $REVERSE_FILE"
echo "  SLURM job name............: $SBATCH_JOB_NAME"
echo "  SLURM CPUs................: $SBATCH_CPUS"
echo "  SLURM memory..............: $SBATCH_MEM"
echo "  SLURM time................: $SBATCH_TIME"
echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
echo "  SLURM mail user...........: $SBATCH_MAILUSER"
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""

# Prepare the changes and reverse options for the Python script
CHANGES_FLAG=""
REVERSE_FLAG=""

if [[ -n "$CHANGES_FILE" ]]; then
    CHANGES_FLAG="--changes $CHANGES_FILE"
fi

if [[ -n "$REVERSE_FILE" ]]; then
    REVERSE_FLAG="--reverse $REVERSE_FILE"
fi

cat << EOF > $SBATCH_SCRIPT
#!/bin/bash
#SBATCH --job-name=$SBATCH_JOB_NAME
#SBATCH --cpus-per-task=$SBATCH_CPUS
#SBATCH --mem=$SBATCH_MEM
#SBATCH --time=$SBATCH_TIME
#SBATCH --mail-type=$SBATCH_MAILTYPE
#SBATCH --mail-user=$SBATCH_MAILUSER
#SBATCH --output=${SBATCH_JOB_NAME}_%j.out
#SBATCH --error=${SBATCH_JOB_NAME}_%j.err

source ~/.bashrc
mamba activate monopogen

echo "$VERSION_NAME"
echo "version $VERSION ($VERSION_DATE)"
echo ""
echo "These are the settings:"
echo "  Input file................: $INPUT_FILE"
echo "  Output file...............: $OUTPUT_FILE"
echo "  Changes file..............: $CHANGES_FILE"
echo "  Reverse file..............: $REVERSE_FILE"
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
python3 $MPG/src/makediploidmalesX.py --input-file $INPUT_FILE --output-file $OUTPUT_FILE ${CHANGES_FLAG} ${REVERSE_FLAG} ${VERBOSE:+--verbose}

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