#!/bin/bash 

# compute_ld_matrices_for_significant_blocks.sh
# Computes LD matrices for significant LD blocks.

# Load necessary modules (adjust as needed)
# module load plink/2.0
# module load plink/1.9
# module load bedtools
# module load R

# Directories and files
TMP_OUTPUT_DIR="tmp_output"
TMP_LD_DIR="tmp_ld"
PGEN_DIR="merged_arrays"  # Update this path
LD_BLOCKS_FILE="ld_blocks/ld_blocks_with_ids.bed"
SIGNIFICANT_LD_BLOCKS_FILE="${TMP_OUTPUT_DIR}/significant_ld_blocks.txt"
OUTPUT_DIR="ld_matrices"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TMP_LD_DIR"

# Read significant block IDs into an array
mapfile -t BLOCK_LINES < "$SIGNIFICANT_LD_BLOCKS_FILE"

# Loop over each significant LD block
for BLOCK_LINE in "${BLOCK_LINES[@]}"; do
    # Extract block information
    CHR=$(echo "$BLOCK_LINE" | awk '{print $1}' | sed 's/chr//')
    START=$(echo "$BLOCK_LINE" | awk '{print $2}')
    END=$(echo "$BLOCK_LINE" | awk '{print $3}')
    BLOCK_ID=$(echo "$BLOCK_LINE" | awk '{print $4}')

    echo "Processing LD block $BLOCK_ID on chromosome $CHR"

    PVAR_FILE="${PGEN_DIR}/gsa_pmra_imputed_chr${CHR}_maf01_rsq8_geno01_ipn_id_excluded_array_assoc.pvar"
#    PVAR_FILE="${PGEN_DIR}/kgp_chr${CHR}.pvar"

    if [ ! -f "$PVAR_FILE" ]; then
        echo "PVAR file for chromosome $CHR not found at $PVAR_FILE. Skipping block $BLOCK_ID."
        continue
    fi

    # Define the expected output LD matrix file
    LD_MATRIX_FILE="${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.ld.gz"
    
    if [ -f "$LD_MATRIX_FILE" ]; then
        echo "LD matrix for block $BLOCK_ID already exists at $LD_MATRIX_FILE. Skipping computation."
        continue  # Skip to the next LD block
    fi


    # Extract SNP positions and IDs from the reference data
    awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1, $2-1, $2, $3}' "$PVAR_FILE" > "${TMP_LD_DIR}/chr${CHR}_snps.bed"

    # Find SNPs in the LD block
    awk -v chr="$CHR" -v start="$START" -v end="$END" \
        'BEGIN{FS=OFS="\t"} ($1 == chr) && ($2 >= start) && ($3 <= end)' "${TMP_LD_DIR}/chr${CHR}_snps.bed" > "${TMP_LD_DIR}/snps_in_block_${BLOCK_ID}.bed"

    # Extract SNP IDs
    cut -f4 "${TMP_LD_DIR}/snps_in_block_${BLOCK_ID}.bed" > "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt"

    # Check if SNPs exist
    if [ ! -s "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt" ]; then
        echo "No SNPs found for block $BLOCK_ID. Skipping."
        rm "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt"
        continue
    fi

    # Extract SNP IDs starting with HGSV
    awk '$3 ~ /^HGSV/ { print $3}' "${PGEN_DIR}/gsa_pmra_imputed_chr${CHR}_maf01_rsq8_geno01_ipn_id_excluded_array_assoc.pvar" > "${TMP_LD_DIR}/exclude_variants_${BLOCK_ID}.txt" \ 

    # Convert PGEN to BED using PLINK 2.0
    plink2 --pfile "${PGEN_DIR}/gsa_pmra_imputed_chr${CHR}_maf01_rsq8_geno01_ipn_id_excluded_array_assoc" \
           --snps-only just-acgt \
           --extract "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt" \
           --exclude "${TMP_LD_DIR}/exclude_variants_${BLOCK_ID}.txt" \
           --make-bed \
           --out "${TMP_LD_DIR}/ref_block_${BLOCK_ID}" \
           --threads 1

# Process the input file
while IFS=$'\t' read -r col1 col2 col3 col4; do
    # Extract values from the second column
    IFS=':' read -r chr pos ref alt <<< "$col2"
    
    # Write the formatted output to the new file
    echo -e "$col2\t$chr\t$pos\t$ref\t$alt" >> "${TMP_LD_DIR}/ref_block_${BLOCK_ID}.txt"
done < "${TMP_LD_DIR}/ref_block_${BLOCK_ID}.bim"

    # Compute LD matrix using PLINK 1.9
    plink --bfile "${TMP_LD_DIR}/ref_block_${BLOCK_ID}" \
          --a1-allele "${TMP_LD_DIR}/ref_block_${BLOCK_ID}.txt" 5 1 '#' \
          --r square \
          --out "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}" \
          --threads 1

    # Save SNP IDs
    cp "${TMP_LD_DIR}/ref_block_${BLOCK_ID}.bim" "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.bim"
    gzip "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.ld"

    # Clean up temporary files for this block
    rm "${TMP_LD_DIR}/chr${CHR}_snps.bed" "${TMP_LD_DIR}/snps_in_block_${BLOCK_ID}.bed" "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt" 
    rm "${TMP_LD_DIR}/ref_block_${BLOCK_ID}".bed "${TMP_LD_DIR}/ref_block_${BLOCK_ID}".bim "${TMP_LD_DIR}/ref_block_${BLOCK_ID}".fam "${TMP_LD_DIR}/ref_block_${BLOCK_ID}".txt "${TMP_LD_DIR}/ref_block_${BLOCK_ID}".log "${TMP_LD_DIR}/exclude_variants_${BLOCK_ID}".txt 

    echo "LD matrix for block $BLOCK_ID computed successfully."
done

