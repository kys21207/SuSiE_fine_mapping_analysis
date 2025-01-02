#!/bin/bash 

# prepare_susie_input_for_significant_blocks.sh
# Prepares input data for SuSiE for significant LD blocks.

# Load necessary modules (adjust as needed)
# module load R

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory and a file name for a GWAS summary statistic."
    echo "Usage: $0 GWAS_SUM_DIR GWAS_SUM_file_name"
    exit 1
fi

# Path to your gzipped GWAS summary statistics file
SUMSTATS_DIR=$1
GWAS_FILE_NAME=$2
SAMPLE_SIZE=$3

# Directories and files
TMP_OUTPUT_DIR="tmp_output"
SUMSTATS_FILE=${SUMSTATS_DIR}/${GWAS_FILE_NAME}
SIGNIFICANT_LD_BLOCKS_FILE="${TMP_OUTPUT_DIR}/significant_ld_blocks.txt"
LD_MATRICES_DIR="ld_matrices"
OUTPUT_DIR="susie_input"
mkdir -p "$OUTPUT_DIR"

# Delimiter: space-delimited
DELIMITER="\t"      # Space delimiter

# Column numbers (adjust if necessary)
CHROM_COLUMN=1    # CHROM
GENPOS_COLUMN=2   # GENPOS
SNPID_COLUMN=5    # ID

# Loop over each significant LD block
while read -r BLOCK_LINE; do
    # Extract block information
    CHR=$(echo "$BLOCK_LINE" | awk '{print $1}' | sed 's/chr//')
    START=$(echo "$BLOCK_LINE" | awk '{print $2}')
    END=$(echo "$BLOCK_LINE" | awk '{print $3}')
    BLOCK_ID=$(echo "$BLOCK_LINE" | awk '{print $4}')

    echo "Preparing SuSiE input for block $BLOCK_ID on chromosome $CHR"

    # Extract summary statistics for SNPs in the block
    gunzip -c "$SUMSTATS_FILE" | awk -v chr_col=$CHROM_COLUMN -v pos_col=$GENPOS_COLUMN -v snp_col=$SNPID_COLUMN \
        -v chr="$CHR" -v start="$START" -v end="$END" -v FS="$DELIMITER" -v OFS="\t" ' NR==1 || ($chr_col == chr && $pos_col >= start && $pos_col <= end)' \
         > "${OUTPUT_DIR}/sumstats_block_${BLOCK_ID}.txt"

    # Prepare data in R
    Rscript codes/prepare_susie_data_cox.R "$BLOCK_ID" "$SAMPLE_SIZE"
done < "$SIGNIFICANT_LD_BLOCKS_FILE"


