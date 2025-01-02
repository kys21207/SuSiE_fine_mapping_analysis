#!/bin/bash 

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a file name for the final results."
    echo "Usage: $0 file_name"
    exit 1
fi

finaloutputname=$1

# run_susie_for_significant_blocks.sh
# Runs SuSiE for significant LD blocks.

TMP_OUTPUT_DIR="tmp_output"
# Load necessary modules (adjust as needed)
# module load R

# Directories and files
SIGNIFICANT_LD_BLOCKS_FILE="${TMP_OUTPUT_DIR}/significant_ld_blocks.txt"

# Loop over each significant LD block
while read -r BLOCK_LINE; do
    # Extract block ID
    BLOCK_ID=$(echo "$BLOCK_LINE" | awk '{print $4}')

    echo "Running SuSiE for block $BLOCK_ID"

    # Run SuSiE in R
    Rscript codes/run_susie_cox.R "$BLOCK_ID" "$finaloutputname"
done < "$SIGNIFICANT_LD_BLOCKS_FILE"

# Clean all SuSiE input files
rm "susie_input/sumstats_block_"*.txt "susie_input/susie_data_block_"*.rds

