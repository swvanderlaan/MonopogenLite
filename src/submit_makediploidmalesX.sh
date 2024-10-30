#!/bin/bash

# Change log:
# * v1.1.8 2024-10-30: Added options to deal with gpu partitions. Added --account option to specify the SLURM account to use. 
# * v1.1.7 2024-10-28: Fixed an issue where a partition can be given when submitting jobs on the UVA RIVANNA cluster.
# * v1.1.6 2024-09-30: Fixed an issue where the concatenated VCF was not done in order.
# * v1.1.5 2024-09-30: Fixed an issue where the concatenated VCF was not sorted prior to indexing.
# * v1.1.4 2024-09-30: Fixed an issue where the chrX was not written correctly to the bed file.
# * v1.1.3 2024-09-30: Fixed an issue where the script was not properly creating the chunks and intermediate variables.
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
VERSION='1.1.8'
VERSION_DATE='2024-10-30'
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
SBATCH_MAILUSER=""
SBATCH_PARTITION="cpu"  # Default partition
# GPU settings
SBATCH_MEM_GPU="16G"
SBATCH_GPUS_NODE="1"
SBATCH_TASKS=1
SBATCH_ACCOUNT="dhl_ec" # Default account
# Other default values
VERBOSE=0  # Set to 0 by default (not verbose)
DEBUG=0  # Set to 0 by default (not debug)
CHUNK_SIZE_DEFAULT=100 # Default number of chunks to process in one go
CHANGES_FILE="" # Default changes file, none
REVERSE_FILE="" # Default reverse file, none
DRY_RUN=0  # Set to 0 by default (not a dry run)

# Argument parsing function
print_help() {
    echo "$VERSION_NAME version $VERSION ($VERSION_DATE)"
    echo ""
    echo "Usage: $0 --input <input.vcf.gz> --output <output.vcf.gz> --mpg-dir [/dir/to/MonopogenLite] [--chunk-size <#>] [--changes <changes.txt.gz>] [--reverse <reverse.txt.gz>] [--job-name <job_name>] [--cpus <num_cpus>] [--tasks <tasks>] [--mem <memory>] [--time <time>] [--partition <type of processor>] [--mem_gpu <memory>] [--gpus_node <num_gpus>] [--account <account>] [--dry-run] [--verbose] [--debug] [--help] [--version]"
    echo ""
    echo "Description:"
    echo "  This script will submit a job to the SLURM scheduler to make haploid genotypes in males diploid given a VCF file."
    echo ""
    echo "Arguments:"
    echo "  --input         The input VCF file. Required."
    echo "  --output        The output VCF file. Required."
    echo "  --mpg-dir       The directory where MonopogenLite is installed. Default is $MPG_DIR. Optional."
    echo "  --chunk-size    The number of chunks to process in one go. Default is 100. Optional."
    echo "  --changes       Gzipped file to save the list of changes. Optional."
    echo "  --reverse       Gzipped file with list of changes to reverse. Optional."
    echo "  --cpus          The number of CPUs to use. Default is 4. Optional."
    echo "  --tasks         The number of tasks to use. Default is 1. Optional."
    echo "  --mem           The amount of memory to use. Default is 16G. Optional."
    echo "  --time          The maximum time to run the job. Default is 1 hour. Optional."
    echo "  --mail          The type of mail to send. Default is FAIL. Optional."
    echo "  --user          The email address to send the mail to. Optional."
    echo "  --partition     SLURM partition to use (default: cpu). Optional."
    echo "  --mem_gpu       The amount of memory to use on the GPU. Default is 16G. Optional."
    echo "  --gpus_node     The number of GPUs per node. Default is 1. Optional."
    echo "  --account       The SLURM account to use. Default is dhl_ec. Optional."
    echo "  --dry-run       Perform a dry run without submitting the job. Optional."
    echo "  --verbose       Enable verbose output. Optional."
    echo "  --debug         Enable debug mode. Optional."
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

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_FILE="$2"; shift ;;
        --output) OUTPUT_FILE="$2"; shift ;;
        --mpg-dir) MPG_DIR="$2"; shift ;;
        --chunk-size) CHUNK_SIZE="$2"; shift ;;
        --changes) CHANGES_FILE="$2"; shift ;;
        --reverse) REVERSE_FILE="$2"; shift ;;
        --cpus) SBATCH_CPUS="$2"; shift ;;
        --tasks) SBATCH_TASKS="$2"; shift ;;
        --mem) SBATCH_MEM="$2"; shift ;;
        --time) SBATCH_TIME="$2"; shift ;;
        --mail) SBATCH_MAILTYPE="$2"; shift ;;
        --user) SBATCH_MAILUSER="$2"; shift ;;
        --partition) PARTITION="$2"; shift ;;
        --mem-gpu) SBATCH_MEM_GPU="$2"; shift ;;
        --gpus-node) SBATCH_GPUS_NODE="$2"; shift ;;
        --account) SBATCH_ACCOUNT="$2"; shift ;;
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
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$MPG_DIR" ]]; then
    echo "Error: The --input, --output --mpg-dir arguments are required."
    echo ""
    print_help
    exit 1
fi

# MonopogenLite location
MPG=$MPG_DIR

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    echo ""
    print_help
    exit 1
fi

# Extract the base name of the input file and its directory
BASE_NAME_INPUT_FILE=$(basename "$INPUT_FILE" .vcf.gz)
BASE_NAME_INPUT_DIR=$(dirname "$INPUT_FILE")
BASE_NAME_OUTPUT_FILE=$(basename "$OUTPUT_FILE" .vcf.gz)

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
SBATCH_SCRIPT_MAKEDIPLOIDMALESX="$MPG/submit_makediploidmalesX.sbatch"
SBATCH_SCRIPT_CONCATINDEX="$MPG/submit_concatindexdiploidmalesX.sbatch"

echo "These are the settings:"
echo "  Input file................: $INPUT_FILE"
echo "  Output file...............: $OUTPUT_FILE"
echo "  MPG directory.............: $MPG_DIR"
echo "  Chunk size................: $CHUNK_SIZE"
if [[ -n "$CHANGES_FILE" ]]; then
    echo "  Changes file..............: $CHANGES_FILE"
fi
if [[ -n "$REVERSE_FILE" ]]; then
    echo "  Reverse file..............: $REVERSE_FILE"
fi
echo "  SLURM CPUs................: $SBATCH_CPUS"
echo "  SLURM tasks...............: $SBATCH_TASKS"
echo "  SLURM partition...........: $PARTITION"
echo "  SLURM time................: $SBATCH_TIME"
if [[ "$PARTITION" == "gpu" ]]; then
    echo "  SLURM mem GPU.............: $SBATCH_MEM_GPU"
    echo "  SLURM GPUs per node.......: $SBATCH_GPUS_NODE"
else
    echo "  SLURM memory..............: $SBATCH_MEM"
fi
echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
echo "  SLURM mail user...........: $SBATCH_MAILUSER"
echo "  SLURM account.............: $SBATCH_ACCOUNT"
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

# Create the BED file with the chunks
for i in $(seq 1 $CHUNK_SIZE); do
    START=$(( (i - 1) * INTERVAL_SIZE + 1 ))
    END=$(( i * INTERVAL_SIZE ))
    if [[ $i -eq $CHUNK_SIZE ]]; then
        END=$CHROM_SIZE  # Ensure last chunk ends at chromosome end
    fi
    echo -e "chrX\t$START\t$END"
done > $BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed

if [[ $DEBUG -eq 1 ]]; then
    echo "DEBUG: BED file created."
    echo "DEBUG: > Check the BED file."
    head $BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed
    echo "DEBUG: > Count the number of lines in the BED file: $(cat $BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed | wc -l)."
fi

echo "> Chunking done."
echo ""

# Create the SLURM batch job script
echo "> Create the SLURM batch job script for chunking haploid genotypes in given VCF file."

# Define SLURM options based on partition type for MakeDiploidMalesX
SBATCH_PARTITION_OPTIONS_MAKEDIPLOID=""
if [[ "$PARTITION" == "gpu" ]]; then
    SBATCH_PARTITION_OPTIONS_MAKEDIPLOID="#SBATCH --gpus-per-node=$SBATCH_GPUS_NODE
#SBATCH --mem-per-gpu=$SBATCH_MEM_GPU
#SBATCH -A $SBATCH_ACCOUNT
#SBATCH --ntasks=$SBATCH_TASKS
#SBATCH --time=$SBATCH_TIME"
else
    SBATCH_PARTITION_OPTIONS_MAKEDIPLOID="#SBATCH --cpus-per-task=$SBATCH_CPUS
#SBATCH --mem=$SBATCH_MEM
#SBATCH --time=$SBATCH_TIME"
fi

cat << EOF > $SBATCH_SCRIPT_MAKEDIPLOIDMALESX
#!/bin/bash
#SBATCH --job-name=makediploidmalesX
#SBATCH --array=1-$CHUNK_SIZE
#SBATCH --partition=$PARTITION
#SBATCH --mail-type=$SBATCH_MAILTYPE
#SBATCH --mail-user=$SBATCH_MAILUSER
#SBATCH --output=makediploidmalesX_%A_%a.out
#SBATCH --error=makediploidmalesX_%A_%a.err
$SBATCH_PARTITION_OPTIONS_MAKEDIPLOID

source ~/.bashrc
source ~/.bash_profile
micromamba activate monopogen

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
echo "  SLURM partition...........: $PARTITION"
echo "  SLURM time................: $SBATCH_TIME"
if [[ $PARTITION == "gpu" ]]; then
    echo "  SLURM mem GPU.............: $SBATCH_MEM_GPU"
    echo "  SLURM GPUs per node.......: $SBATCH_GPUS_NODE"
    echo "  SLURM tasks...............: $SBATCH_TASKS"
    echo "  SLURM account.............: $SBATCH_ACCOUNT"
else 
    echo "  SLURM memory..............: $SBATCH_MEM"
fi
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

echo "  - Get the chunk number from the SLURM array task ID and the bed file..."
CHUNK_NUMBER_RAW=\$SLURM_ARRAY_TASK_ID # SLURM array task ID
CHUNK_NUMBER=\${CHUNK_NUMBER_RAW}p  # adding p to the number
BED_FILE="$BASE_NAME_INPUT_DIR/chrX_${CHUNK_SIZE}pieces.bed" # BED file with chunks

# Debugs
if [[ \$DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Extracting chunk number: \$CHUNK_NUMBER from \$CHUNK_NUMBER_RAW"
fi
if [[ \$DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Extracting region from BED file: \$BED_FILE"
fi

echo "  - Extract the region from the BED file..."
RAW_REGION=\$(sed -n "\$CHUNK_NUMBER" "\$BED_FILE")
RAW_REGION=\$(echo "\$RAW_REGION" | sed 's/^ *//;s/ *\$//')  # Clean up leading/trailing whitespace
if [[ \$DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Extracted raw region: '\$RAW_REGION' (chunk #: \$CHUNK_NUMBER of $CHUNK_SIZE for job \$CHUNK_NUMBER_RAW)"
fi
echo "  - Check if a valid region was extracted -- should be '\$RAW_REGION'..."
if [[ -z "\$RAW_REGION" ]]; then
  echo "ERROR: No valid region extracted for chunk # \$SLURM_ARRAY_TASK_ID (or \$CHUNK_NUMBER_RAW) from '\$BED_FILE'."
  exit 1
fi

echo "  - Prepare the region for bcftools..."
REGION=\$(echo "\$RAW_REGION" | awk '{print \$1 ":" \$2 "-" \$3}')
echo "  - Check if the region (\$REGION) was successfully extracted..."
if [[ -z "\$REGION" ]]; then
  echo "ERROR: RAW_REGION was: '\$RAW_REGION'"
  echo "ERROR: No region found for SLURM_ARRAY_TASK_ID \$SLURM_ARRAY_TASK_ID."
  exit 1
fi

echo "> Processing region '\$REGION'..."

echo "  - Process the specific region (chunk# \$CHUNK_NUMBER_RAW) using bcftools view..."
bcftools view -r \$REGION $INPUT_FILE -Oz -o $BASE_NAME_INPUT_DIR/chrX.part\${CHUNK_NUMBER_RAW}.vcf.gz

echo "  - Index the VCF file..."
tabix -fp vcf $BASE_NAME_INPUT_DIR/chrX.part\${CHUNK_NUMBER_RAW}.vcf.gz

# Debugs
if [[ \$DEBUG_FLAG -eq "$DEBUG" ]]; then
    echo "DEBUG: Check if the VCF file was successfully processed..."
    bcftools view $BASE_NAME_INPUT_DIR/chrX.part\${CHUNK_NUMBER_RAW}.vcf.gz | head
    bcftools stats $BASE_NAME_INPUT_DIR/chrX.part\${CHUNK_NUMBER_RAW}.vcf.gz
fi

echo "  - Fix haploid genotypes in males..."
python3 $MPG/src/makediploidmalesX.py --input-file $BASE_NAME_INPUT_DIR/chrX.part\${CHUNK_NUMBER_RAW}.vcf.gz --output-file $BASE_NAME_INPUT_DIR/chrX.part\${CHUNK_NUMBER_RAW}_processed.vcf.gz ${CHANGES_FLAG} ${REVERSE_FLAG} ${VERBOSE_FLAG}

if [ \$? -eq 0 ]; then
    echo ""
    echo "$VERSION_NAME finished successfully. Let's have a beer, buddy!"
    echo ""
else
    echo ""
    echo "ERROR: $VERSION_NAME encountered an error."
    echo ""
fi

micromamba deactivate

EOF

# Make the script executable
echo "> Make the script executable."
chmod +x $SBATCH_SCRIPT_MAKEDIPLOIDMALESX

# Create the SLURM batch job script
echo "> Create the SLURM batch job script for concatenating the processed VCF files."

# Define SLURM options based on partition type for ConcatIndex job
SBATCH_PARTITION_OPTIONS_CONCAT=""
if [[ "$PARTITION" == "gpu" ]]; then
    SBATCH_PARTITION_OPTIONS_CONCAT="#SBATCH --gpus-per-node=$SBATCH_GPUS_NODE
#SBATCH --mem-per-gpu=128G
#SBATCH -A $SBATCH_ACCOUNT
#SBATCH --ntasks=$SBATCH_TASKS
#SBATCH --time=02:00:00"
else
    SBATCH_PARTITION_OPTIONS_CONCAT="#SBATCH --cpus-per-task=$SBATCH_CPUS
#SBATCH --mem=128G
#SBATCH --time=02:00:00"
fi


cat << EOF > $SBATCH_SCRIPT_CONCATINDEX
#!/bin/bash
#SBATCH --job-name=concatdiploidmalesX
#SBATCH --partition=$PARTITION
#SBATCH --mail-type=$SBATCH_MAILTYPE
#SBATCH --mail-user=$SBATCH_MAILUSER
#SBATCH --output=concatdiploidmalesX_%j.out
#SBATCH --error=concatdiploidmalesX_%j.err
$SBATCH_PARTITION_OPTIONS_CONCAT

source ~/.bashrc
source ~/.bash_profile
micromamba activate monopogen

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
echo "  SLURM partition...........: $PARTITION"
echo "  SLURM time................: $SBATCH_TIME"
if [[ $PARTITION == "gpu" ]]; then
    echo "  SLURM mem GPU.............: $SBATCH_MEM_GPU"
    echo "  SLURM GPUs per node.......: $SBATCH_GPUS_NODE"
    echo "  SLURM tasks...............: $SBATCH_TASKS"
    echo "  SLURM account.............: $SBATCH_ACCOUNT"
else 
    echo "  SLURM memory..............: $SBATCH_MEM"
fi
echo "  SLURM mail type...........: $SBATCH_MAILTYPE"
echo "  SLURM mail user...........: $SBATCH_MAILUSER"
echo ""
echo "  Dry run mode..............: $DRY_RUN"
echo "  Debug mode................: $DEBUG"
echo "  Verbosity.................: $VERBOSE"
echo "  Version...................: $VERSION ($VERSION_DATE)"
echo ""
echo "Concatenating processed VCF files."

echo "> Concatenate the processed VCF files..."
echo "  - Construct the ordered list of files."
VCF_LIST=""
CHUNKSIZE=$CHUNK_SIZE
for CHUNK in \$(seq 1 \$CHUNKSIZE); do
    FILE="$BASE_NAME_INPUT_DIR/chrX.part\${CHUNK}_processed.vcf.gz"
    echo "  - Processing chunk \$CHUNK and adding [\$FILE] to the list..."
    if [[ -f "\$FILE" ]]; then
        VCF_LIST="\$VCF_LIST \$FILE"
    else
        echo "ERROR: File \$FILE does not exist..."
        exit 1
    fi
done

if [[ \$DEBUG_FLAG -eq 1 ]]; then
    echo "DEBUG: Listed all the processed VCF files:"
    printf "%s\n" "\$VCF_LIST"
fi

bcftools concat -Oz -o $OUTPUT_FILE \$VCF_LIST
if [ $? -ne 0 ]; then
    echo "ERROR: Concatenation failed."
    exit 1
fi

echo "> Index the concatenated VCF file..."
tabix -fp vcf $OUTPUT_FILE
if [ $? -ne 0 ]; then
    echo "ERROR: Indexing with tabix failed."
    exit 1
fi

echo "> Removing the processed VCF files..."
for CHUNK in \$(seq 1 \$CHUNKSIZE); do
    FILE="$BASE_NAME_INPUT_DIR/chrX.part\${CHUNK}_processed.vcf.gz"
    FILE_TBI="$BASE_NAME_INPUT_DIR/chrX.part\${CHUNK}_processed.vcf.gz.tbi"
    FILE_RAW="$BASE_NAME_INPUT_DIR/chrX.part\${CHUNK}.vcf.gz"
    FILE_RAW_TBI="$BASE_NAME_INPUT_DIR/chrX.part\${CHUNK}.vcf.gz.tbi"
    echo "  - Removing processed VCF file: [\$FILE]"
    rm -v \$FILE
    echo "  - Removing processed VCF file: [\$FILE_TBI]"
    rm -v \$FILE_TBI
    echo "  - Removing raw VCF file: [\$FILE_RAW]"
    rm -v \$FILE_RAW
    echo "  - Removing raw VCF file: [\$FILE_RAW_TBI]"
    rm -v \$FILE_RAW_TBI
done

if [ \$? -eq 0 ]; then
    echo ""
    echo "$VERSION_NAME finished successfully. Let's have a beer, buddy!"
    echo ""
else
    echo ""
    echo "ERROR: $VERSION_NAME encountered an error."
    echo ""
fi

micromamba deactivate

EOF

# Make the script executable
echo "> Make the script executable."
chmod +x $SBATCH_SCRIPT_CONCATINDEX

echo ""
echo "> Submit the array job, dependent on the BED creation job."
if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN: sbatch $SBATCH_SCRIPT_MAKEDIPLOIDMALESX"
else
    ARRAY_JOB_ID=$(sbatch $SBATCH_SCRIPT_MAKEDIPLOIDMALESX | awk '{print $4}')
    echo ">> Array job submitted with ID: $ARRAY_JOB_ID."
fi

echo "> Submit the concatenation job, dependent on the array job completion."
if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN: sbatch --dependency=afterok:$ARRAY_JOB_ID $SBATCH_SCRIPT_CONCATINDEX"
else
    CONCAT_JOB_ID=$(sbatch --dependency=afterok:$ARRAY_JOB_ID $SBATCH_SCRIPT_CONCATINDEX | awk '{print $4}')
    echo ">> Concatenation job submitted with ID: $CONCAT_JOB_ID."
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
