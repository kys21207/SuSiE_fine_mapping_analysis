#!/bin/bash 

# filter_sumstats.sh
# Filters GWAS summary statistics for SNPs with LOG10P > 5 (p-value < 1e-5).
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

LD_BLOCKS_FILE="ld_blocks/ld_blocks_with_ids.bed"
GWAS_ID=$(basename "$GWAS_FILE_NAME" .tsv.gz)
SUMSTATS_FILE=${SUMSTATS_DIR}/${GWAS_FILE_NAME}
TMP_DIR="tmp"
FILTERED_SUMSTATS_FILE="${TMP_DIR}/${GWAS_ID}_sumstats_pvalue_lt_1e8.txt"
# Output file
SIGNIFICANT_LD_BLOCKS_FILE="${TMP_DIR}/${GWAS_ID}_significant_ld_blocks.txt"
SUSIE_INPUT_DIR="susie_input/${GWAS_ID}"

mkdir -p "$TMP_DIR"
mkdir -p "$SUSIE_INPUT_DIR"
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE CHISQ LOG10P INFO

# Column number for LOG10P (1-based indexing)
LOG10P_COLUMN=11   # Adjust if LOG10P is in a different column

# Delimiter: space-delimited
DELIMITER=" "      # Space delimiter
#DELIMITER="\t"      # Tab delimiter

# Check if the summary statistics file exists
if [ ! -f "$SUMSTATS_FILE" ]; then
    echo "Error: Summary statistics file '$SUMSTATS_FILE' not found."
    exit 1
fi

# Determine if the file is space-delimited by inspecting the first data line
FIRST_LINE=$(gunzip -c "$SUMSTATS_FILE" | grep -v "^CHROM" | head -n 1)
#FIRST_LINE=$(tail -n +2 "$SUMSTATS_FILE" | grep -v "^CHROM" | head -n 1)

# Count the number of fields in the first data line
NUM_FIELDS=$(echo "$FIRST_LINE" | awk '{print NF}')

echo "Detected $NUM_FIELDS fields in the summary statistics file."

# Verify that LOG10P_COLUMN does not exceed the number of fields
if [ "$LOG10P_COLUMN" -gt "$NUM_FIELDS" ]; then
    echo "Error: LOG10P_COLUMN ($LOG10P_COLUMN) exceeds the number of fields ($NUM_FIELDS)."
    exit 1
fi

# Filter SNPs with LOG10P > 5 and ensure LOG10P is numeric
echo "Filtering SNPs with LOG10P > 8..."

gunzip -c "$SUMSTATS_FILE" | awk -v log10p_col=$LOG10P_COLUMN -v FS="$DELIMITER" -v OFS="\t" '
NR == 1 {print; next} 
{
    # Check if LOG10P is a number and greater than 8
    if ($log10p_col ~ /^-?[0-9]*\.?[0-9]+$/ && $log10p_col > 8) {
        print
    }
}' > "$FILTERED_SUMSTATS_FILE"

# Count the number of SNPs after filtering
TOTAL_SNPS=$(gunzip -c "$SUMSTATS_FILE" | wc -l)
#TOTAL_SNPS=$(wc -l < "$SUMSTATS_FILE")
FILTERED_SNPS=$(wc -l < "$FILTERED_SUMSTATS_FILE")

echo "Total SNPs in original file: $TOTAL_SNPS"
echo "Total SNPs after filtering: $FILTERED_SNPS"

# Optional: Report how many SNPs were filtered out
SNPS_FILTERED=$((TOTAL_SNPS - FILTERED_SNPS))
echo "Number of SNPs filtered out: $SNPS_FILTERED"

# Column numbers (adjust if necessary)
CHROM_COLUMN=1    # CHROM is column 1
GENPOS_COLUMN=2   # GENPOS is column 2
SNPID_COLUMN=3    # ID is column 3

# Convert summary statistics to BED format
awk -v chr_col=$CHROM_COLUMN -v pos_col=$GENPOS_COLUMN -v snp_col=$SNPID_COLUMN -v FS="$DELIMITER" -v OFS="\t" \
    'NR>1 {print "chr"$chr_col, $pos_col-1, $pos_col, $snp_col}' "$FILTERED_SUMSTATS_FILE" > "${TMP_DIR}/${GWAS_ID}_significant_snps.bed"

bedtools intersect -a "${LD_BLOCKS_FILE}" -b "${TMP_DIR}/${GWAS_ID}_significant_snps.bed" -wa | sort | uniq > "${SIGNIFICANT_LD_BLOCKS_FILE}"



