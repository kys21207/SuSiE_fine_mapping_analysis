#!/bin/bash 
# Computes LD matrices for all LD blocks.
# Directories and files
TMP_OUTPUT_DIR="tmp_output"
TMP_LD_DIR="tmp_ld"
PGEN_DIR="/mnt/project/10k_ld_reference_panel_topmed"  # Update this path
LD_BLOCKS_FILE="ld_blocks/ld_blocks_with_ids.bed"
#SIGNIFICANT_LD_BLOCKS_FILE="${TMP_OUTPUT_DIR}/significant_ld_blocks.txt"
OUTPUT_DIR="ld_matrices"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TMP_LD_DIR"

# Read significant block IDs into an array
mapfile -t BLOCK_LINES < "$LD_BLOCKS_FILE"

# Loop over each significant LD block
for BLOCK_LINE in "${BLOCK_LINES[@]}"; do
    # Extract block information
    CHR=$(echo "$BLOCK_LINE" | awk '{print $1}' | sed 's/chr//')
    START=$(echo "$BLOCK_LINE" | awk '{print $2}')
    END=$(echo "$BLOCK_LINE" | awk '{print $3}')
    BLOCK_ID=$(echo "$BLOCK_LINE" | awk '{print $4}')

    # Skip blocks outside the target range
    if [ "$BLOCK_ID" -lt 224 ] || [ "$BLOCK_ID" -gt 446 ]; then
        continue
    fi
    
    echo "Processing LD block $BLOCK_ID on chromosome $CHR"

    PVAR_FILE="${PGEN_DIR}/ukb21007_c${CHR}_b0_v1_imputed_maf001_rsq7_10k.pvar"
#    PVAR_FILE="${PGEN_DIR}/kgp_chr${CHR}.pvar"
    FILE_PREFIX=$(basename "$PVAR_FILE" .pvar)

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
    
    # Convert PGEN to BED using PLINK 2.0
    plink2 --pfile "${PGEN_DIR}/${FILE_PREFIX}" \
           --snps-only just-acgt \
           --extract "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt" \
           --r2-unphased square \
           --out "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}" \
           --threads 1

    # Save SNP IDs
    mv "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.vcor2.vars" "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.ld"
    gzip "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.ld"

    dx upload ${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.* --destination project-Gz0KVk8JPx8q9yBzpgXX6y8K:/analysis_supporting_data/ld_matrices/

    # Clean up temporary files for this block
    rm "${TMP_LD_DIR}/snps_in_block_${BLOCK_ID}.bed" "${TMP_LD_DIR}/snps_block_${BLOCK_ID}.txt" "${OUTPUT_DIR}/ld_matrix_block_${BLOCK_ID}.vars"

    echo "LD matrix for block $BLOCK_ID computed successfully."
done
