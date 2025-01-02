#!/bin/bash 

# filter_sumstats.sh
# Filters GWAS summary statistics for SNPs with P-value > 1e-5 (p-value < 1e-5).
# Handles input files compressed in .gz format and space-delimited.

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory and a file name for a GWAS summary statistic."
    echo "Usage: $0 GWAS_SUM_DIR GWAS_SUM_file_name"
    exit 1
fi

# Input and output files
SUMSTATS_DIR=$1
GWAS_FILE_NAME=$2            # Path to your gzipped GWAS summary statistics file
FILTERED_SUMSTATS_FILE="tmp_input/sumstats_pvalue_lt_1e5.txt"
SUMSTATS_FILE=${SUMSTATS_DIR}/${GWAS_FILE_NAME}

# Column number for P-value (1-based indexing)
P_COLUMN=8   # Adjust if P-value is in a different column

# Delimiter: space-delimited
DELIMITER="\t"      # Space delimiter

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$FILTERED_SUMSTATS_FILE")"

# Check if the summary statistics file exists
if [ ! -f "$SUMSTATS_FILE" ]; then
    echo "Error: Summary statistics file '$SUMSTATS_FILE' not found."
    exit 1
fi

# CHROM	GENPOS	ALLELE0	ALLELE1	ID	MAF	missing.rate	p.value.spa	p.value.norm	Stat	Var	z

# Determine if the file is space-delimited by inspecting the first data line
FIRST_LINE=$(gunzip -c "$SUMSTATS_FILE" | grep -v "^CHROM" | head -n 1)

# Count the number of fields in the first data line
NUM_FIELDS=$(echo "$FIRST_LINE" | awk '{print NF}')

echo "Detected $NUM_FIELDS fields in the summary statistics file."

# Verify that LOG10P_COLUMN does not exceed the number of fields
if [ "$P_COLUMN" -gt "$NUM_FIELDS" ]; then
    echo "Error: Pvalue_COLUMN ($P_COLUMN) exceeds the number of fields ($NUM_FIELDS)."
    exit 1
fi

# Filter SNPs with P-value > 1e-5 and ensure P-value is numeric
echo "Filtering SNPs with P-value < 1e-5..."

gunzip -c "$SUMSTATS_FILE" | awk -v pvalue_col="$P_COLUMN" -v FS="$DELIMITER" -v OFS="\t" '
NR == 1 {print; next} 
{
    # Check if P-value is a number and less than 1e-5
    if ($pvalue_col ~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && $pvalue_col < 1e-5) {
        print
    }
}' > "$FILTERED_SUMSTATS_FILE"

# Count the number of SNPs after filtering
TOTAL_SNPS=$(gunzip -c "$SUMSTATS_FILE" | wc -l)
FILTERED_SNPS=$(wc -l < "$FILTERED_SUMSTATS_FILE")

echo "Total SNPs in original file: $TOTAL_SNPS"
echo "Total SNPs after filtering: $FILTERED_SNPS"

# Optional: Report how many SNPs were filtered out
SNPS_FILTERED=$((TOTAL_SNPS - FILTERED_SNPS))
echo "Number of SNPs filtered out: $SNPS_FILTERED"


