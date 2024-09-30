#!/bin/bash

# Change log:
# * v1.1.2 2024-09-30: Added --debug mode. 
# * v1.1.1 2024-09-30: Fixed an issue where the chunking and processing jobs were not properly linked. Added a --dry-run argument to test the script without submitting jobs.
# * v1.1.0 2024-09-27: Added a --chunk-size argument to specify the number of chunks to process in one go thus speeding up the process.
# * v1.0.6 2024-09-27: Fixed an issue where the script was always using --verbose, even if it was not passed.
# * v1.0.5 2024-09-27: Changed default values for SLURM.
# * v1.0.4 2024-09-25: Added a check if the input file exists. Added optional --changes and --reverse flags.
# * v1.0.3 2024-09-25: Changed script name.
# * v1.0.2 2024-09-24: Fixed issue where there was no --version flag in the help message.
# * v1.0.1 2024-09-24: Added a filter for variants that are homozygous in the phased high-coverage 1000 Genomes VCF files, only for chromosome X.
# * v1.0.0 2024-09-24: Initial version. 
# Version and license information 
VERSION_NAME='Submit MakeDiploidMalesX'
VERSION='1.1.2'
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

echo "Starting $VERSION_NAME"
echo ""

# Default values
SBATCH_CPUS=4
SBATCH_MEM="16G"
SBATCH_TIME="01:00:00"
SBATCH_MAILTYPE="FAIL"
SBATCH_MAILUSER="s.w.vanderlaan-2@umcutrecht.nl"
VERBOSE=0  # Set to 0 by default (not verbose)
DEBUG=0  # Set to 0 by default (not debug)
CHUNK_SIZE_DEFAULT=100 # Default number of chunks to process in one go
CHANGES_FILE="" # Default changes file, none
REVERSE_FILE="" # Default reverse file, none
DRY_RUN=0  # Set to 0 by default (not a dry run)

# MonopogenLite location
MPG="/hpc/local/Rocky8/dhl_ec/software/MonopogenLite"

# Argument parsing function
print_help() {
    echo "$VERSION_NAME version $VERSION ($VERSION_DATE)"
    echo ""
    echo "Usage: $0 --input <input.vcf.gz> --output <output.vcf.gz> [--chunk-size <#>] [--changes <changes.txt.gz>] [--reverse <reverse.txt.gz>] [--job-name <job_name>] [--cpus <num_cpus>] [--mem <memory>] [--time <time>] [--mail <mail-type>] [--user <mail-user>] [--verbose]"
    echo ""
    echo "Description:"
    echo "  This script will submit a job to the SLURM scheduler to make haploid genotypes in males diploid given a VCF file."
    echo ""
    echo "Arguments:"
    echo "  --input         The input VCF file."
    echo "  --output        The output VCF file."
    echo "  --chunk-size    The number of chunks to process in one go. Default is 100. Optional."
    echo "  --changes       Gzipped file to save the list of changes (optional)."
    echo "  --reverse       Gzipped file with list of changes to reverse (optional)."
    echo "  --cpus          The number of CPUs to use."
    echo "  --mem           The amount of memory to use."
    echo "  --time          The maximum time to run the job."
    echo "  --mail          The type of mail to send."
    echo "  --user          The email address to send the mail to."
    echo "  --dry-run       Perform a dry run without submitting the job."
    echo "  --verbose       Enable verbose output."
    echo "  --debug         Enable debug mode."
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

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_FILE="$2"; shift ;;
        --output) OUTPUT_FILE="$2"; shift ;;
        --chunk-size) CHUNK_SIZE="$2"; shift ;;
        --changes) CHANGES_FILE="$2"; shift ;;
        --reverse) REVERSE_FILE="$2"; shift ;;
        --cpus) SBATCH_CPUS="$2"; shift ;;
        --mem) SBATCH_MEM="$2"; shift ;;
        --time) SBATCH_TIME="$2"; shift ;;
        --mail) SBATCH_MAILTYPE="$2"; shift ;;
        --user) SBATCH_MAILUSER="$2"; shift ;;
        --dry-run) DRY_RUN=1 ;;  # Set dry run to 1 if --dry-run is passed
        --verbose) VERBOSE=1 ;;  # Set verbose to 1 if --verbose is passed
        --debug) DEBUG=1 ;;  # Set debug to 1 if --debug is passed
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

# Extract the base name of the input file and its directory
BASE_NAME_INPUT_FILE=$(basename "$INPUT_FILE" .vcf.gz)
BASE_NAME_INPUT_DIR=$(dirname "$INPUT_FILE")

# Ensure the input-file directory exists
if [ ! -d "$BASE_NAME_INPUT_DIR" ]; then
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Creating input file directory: $BASE_NAME_INPUT_DIR"
    fi
    mkdir -vp "$BASE_NAME_INPUT_DIR"
fi
# Ensure the output-file directory exists
if [ ! -d "$(dirname $OUTPUT_FILE)" ]; then
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Creating output file directory: $(dirname $OUTPUT_FILE)"
    fi
    mkdir -vp "$(dirname $OUTPUT_FILE)"
fi

# Prepare the chunk size -- set to 100 if not provided or smaller than 2
if [[ -n "$CHUNK_SIZE" && "$CHUNK_SIZE" -ge 2 ]]; then
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Chunk size is set to $CHUNK_SIZE."
    fi
    CHUNK_SIZE=$CHUNK_SIZE
else
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Chunk size is set to default ($CHUNK_SIZE_DEFAULT)."
    fi
    CHUNK_SIZE=$CHUNK_SIZE_DEFAULT
fi

# Prepare the changes and reverse options for the Python script
if [[ -n "$CHANGES_FILE" ]]; then
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Changes file is set to $CHANGES_FILE."
    fi
    CHANGES_FLAG="--changes $CHANGES_FILE"
fi

if [[ -n "$REVERSE_FILE" ]]; then
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Reverse file is set to $REVERSE_FILE."
    fi
    REVERSE_FLAG="--reverse $REVERSE_FILE"
fi

# Only add --verbose if VERBOSE is 1
VERBOSE_FLAG=""
if [[ "$VERBOSE" -eq 1 ]]; then
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "DEBUG: Verbose mode is enabled."
    fi
    VERBOSE_FLAG="--verbose"
fi

# Create the SLURM batch job script
SBATCH_SCRIPT_CHUNKHAPLOIDMALESX="$MPG/submit_chunkhaploidmalesX.sbatch"
SBATCH_SCRIPT_MAKEDIPLOIDMALESX="$MPG/submit_makediploidmalesX.sbatch"
SBATCH_SCRIPT_CONCATINDEX="$MPG/submit_concatindexdiploidmalesX.sbatch"

echo "These are the settings:"
echo "  Input file................: $INPUT_FILE"
echo "  Output file...............: $OUTPUT_FILE"
echo "  Chunk size................: $CHUNK_SIZE"
if [[ -n "$CHANGES_FILE" ]]; then
    echo "  Changes file..............: $CHANGES_FILE"
fi
if [[ -n "$REVERSE_FILE" ]]; then
    echo "  Reverse file..............: $REVERSE_FILE"
fi
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

echo ""
echo "Creating the necessary number of chunks ($CHUNK_SIZE) to split."

echo "> Set chromosome X size (default: 156040895 in GRCh38)."
CHROM_SIZE=156040895  # Length of chromosome X in GRCh38
CHUNK_SIZE_NUMBER=$CHUNK_SIZE        # Number of chunks to split into

INTERVAL_SIZE=$((CHROM_SIZE / $CHUNK_SIZE))

for i in $(seq 1 $CHUNK_SIZE); do
    START=$(( (i - 1) * INTERVAL_SIZE + 1 ))
    END=$(( i * INTERVAL_SIZE ))
    if [[ $i -eq $CHUNK_SIZE ]]; then
        END=$CHROM_SIZE  # Ensure last chunk ends at chromosome end
    fi
    echo -e "X\t$START\t$END"
done > $BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed

echo "> Check the BED file."
head $BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed
cat $BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed | wc -l

echo "Chunking done. Let's process the chunks. This may take a while. Grab a coffee."

# Create the SLURM batch job script
cat << EOF > $SBATCH_SCRIPT_MAKEDIPLOIDMALESX
#!/bin/bash
#SBATCH --job-name=makediploidmalesX
#SBATCH --array=1-$CHUNK_SIZE
#SBATCH --cpus-per-task=$SBATCH_CPUS
#SBATCH --mem=$SBATCH_MEM
#SBATCH --time=$SBATCH_TIME
#SBATCH --mail-type=$SBATCH_MAILTYPE
#SBATCH --mail-user=$SBATCH_MAILUSER
#SBATCH --output=makediploidmalesX_%A_%a.out
#SBATCH --error=makediploidmalesX_%A_%a.err

source ~/.bashrc
mamba activate monopogen

DEBUG_FLAG=$DEBUG

echo "$VERSION_NAME"
echo "version $VERSION ($VERSION_DATE)"
echo ""
echo "These are the settings:"
echo "  Input file................: $INPUT_FILE"
echo "  Output file...............: $OUTPUT_FILE"
echo "  Chunk size................: $CHUNK_SIZE"
if [[ -n "$CHANGES_FILE" ]]; then
    echo "  Changes file..............: $CHANGES_FILE"
fi
if [[ -n "$REVERSE_FILE" ]]; then
    echo "  Reverse file..............: $REVERSE_FILE"
fi
echo ""
echo "  SLURM CPUs................: $SBATCH_CPUS"
echo "  SLURM Array ID............: \$SLURM_ARRAY_TASK_ID (of $CHUNK_SIZE)"
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
echo "Making converting haploid genotypes to diploid genotypes."

echo "> Extract the region for this SLURM task (# \$SLURM_ARRAY_TASK_ID) from the BED file."

# Get the chunk number from the SLURM array task ID and the bed file
CHUNK_NUMBER_RAW=\$SLURM_ARRAY_TASK_ID # SLURM array task ID
CHUNK_NUMBER=\${CHUNK_NUMBER_RAW}p  # adding p to the number
BED_FILE="$BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed" # BED file with chunks

# Debugs
if [[ \$DEBUG_FLAG -eq "$DEBUG" ]]; then
    echo "DEBUG: Extracting chunk number: \$CHUNK_NUMBER from \$CHUNK_NUMBER_RAW"
fi
if [[ \$DEBUG_FLAG -eq "$DEBUG" ]]; then
    echo "DEBUG: Extracting region from BED file: \$BED_FILE"
fi

# Extract the region from the BED file
RAW_REGION=$(sed -n "\$CHUNK_NUMBER" "\$BED_FILE")
if [[ \$DEBUG_FLAG -eq "$DEBUG" ]]; then
    echo "DEBUG: Extracted raw region: \$RAW_REGION (chunk number: \$CHUNK_NUMBER of $CHUNK_SIZE)"
fi
# Prepare the region for bcftools
REGION=$(echo \$RAW_REGION | awk '{print $1 ":" $2 "-" $3}')

if [[ -z \$REGION ]]; then
  echo "Error: No region found for SLURM_ARRAY_TASK_ID \$SLURM_ARRAY_TASK_ID."
  exit 1
fi

echo "> Processing region $REGION..."

echo "  - Process the specific region (chunk) using bcftools view..."
bcftools view -r $REGION $INPUT_FILE -Oz -o $BASE_NAME_INPUT_DIR/chrX.part${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "  - Index the VCF file..."
tabix -p vcf $BASE_NAME_INPUT_DIR/chrX.part${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "  - Fix haploid genotypes in males..."
python3 $MPG/src/makediploidmalesX.py --input-file $BASE_NAME_INPUT_DIR/chrX.part${SLURM_ARRAY_TASK_ID}.vcf.gz --output-file $BASE_NAME_INPUT_DIR/chrX.part${SLURM_ARRAY_TASK_ID}_processed.vcf.gz ${CHANGES_FLAG} ${REVERSE_FLAG} ${VERBOSE_FLAG}

if [ \$? -eq 0 ]; then
  echo "$VERSION_NAME finished successfully. Let's have a beer, buddy!"
else
  echo "$VERSION_NAME encountered an error."
fi

mamba deactivate

EOF

# Make the script executable
chmod +x $SBATCH_SCRIPT_MAKEDIPLOIDMALESX

# # Create the SLURM batch job script
# cat << EOF > $SBATCH_SCRIPT_CONCATINDEX
# #!/bin/bash
# #SBATCH --job-name=concatdiploidmalesX
# #SBATCH --cpus-per-task=$SBATCH_CPUS
# #SBATCH --mem=$SBATCH_MEM
# #SBATCH --time=$SBATCH_TIME
# #SBATCH --mail-type=$SBATCH_MAILTYPE
# #SBATCH --mail-user=$SBATCH_MAILUSER
# #SBATCH --output=concatdiploidmalesX_%j.out
# #SBATCH --error=concatdiploidmalesX_%j.err

# source ~/.bashrc
# mamba activate monopogen

# echo "$VERSION_NAME"
# echo "version $VERSION ($VERSION_DATE)"
# echo ""
# echo "These are the settings:"
# echo "  Input file................: $INPUT_FILE"
# echo "  Output file...............: $OUTPUT_FILE"
# echo "  Chunk size................: $CHUNK_SIZE"
# if [[ -n "$CHANGES_FILE" ]]; then
#     echo "  Changes file..............: $CHANGES_FILE"
# fi
# if [[ -n "$REVERSE_FILE" ]]; then
#     echo "  Reverse file..............: $REVERSE_FILE"
# fi
# echo ""
# echo "  SLURM CPUs................: $SBATCH_CPUS"
# echo "  SLURM memory..............: $SBATCH_MEM"
# echo "  SLURM time................: $SBATCH_TIME"
# echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
# echo "  SLURM mail user...........: $SBATCH_MAILUSER"
# echo ""
# echo "  Dry run mode..............: $DRY_RUN"
# echo "  Debug mode................: $DEBUG"
# echo "  Verbosity.................: $VERBOSE"
# echo "  Version...................: $VERSION ($VERSION_DATE)"
# echo ""
# echo "Concatenating processed VCF files."

# echo "> Concatenate the processed VCF files..."
# bcftools concat -Oz -o $OUTPUT_FILE $BASE_NAME_INPUT_DIR/chrX.part*_processed.vcf.gz

# echo "> Index the concatenated VCF file..."
# bcftools index $OUTPUT_FILE

# if [ \$? -eq 0 ]; then
#   echo "$VERSION_NAME finished successfully. Let's have a beer, buddy!"
# else
#   echo "$VERSION_NAME encountered an error."
# fi

# mamba deactivate

# EOF

# # Make the script executable
# chmod +x $SBATCH_SCRIPT_CONCATINDEX

# Submit the array job, dependent on the BED creation job
if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN: sbatch $SBATCH_SCRIPT_MAKEDIPLOIDMALESX"
else
    ARRAY_JOB_ID=$(sbatch $SBATCH_SCRIPT_MAKEDIPLOIDMALESX | awk '{print $4}')
    echo ">> Array job submitted with ID: $ARRAY_JOB_ID."
fi


# # Submit the concatenation job, dependent on the array job completion
# if [[ $DRY_RUN -eq 1 ]]; then
#     echo "DRY-RUN: sbatch --dependency=afterok:$ARRAY_JOB_ID $SBATCH_SCRIPT_CONCATINDEX"
# else
#     CONCAT_JOB_ID=$(sbatch --dependency=afterok:$ARRAY_JOB_ID $SBATCH_SCRIPT_CONCATINDEX | awk '{print $4}')
#     echo ">> Concatenation job submitted with ID: $CONCAT_JOB_ID."
# fi

echo ""
print_version
### END OF SCRIPT ###
