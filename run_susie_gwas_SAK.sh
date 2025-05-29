#!/bin/bash 

# Directories and files

GWAS_DATA_PATH="/mnt/project/analyses_KJ/meta_EUR_only_results/meta/regenie_format"
OUTPUT_DIR="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/susie_finemapping_analysis/meta_traits"

#for GWAS_FILE in "$GWAS_DATA_PATH"/*.tsv.gz; do
GWAS_FILE="/mnt/project/analyses_KJ/meta_EUR_only_results/meta/regenie_format/G6_MS_meta_out_EUR.tsv.gz"
    GWAS_ID=$(basename "$GWAS_FILE" _meta_out_EUR.tsv.gz)
    GWAS_NAME=$(basename "$GWAS_FILE" )

    # Prepare data in R
    echo "Processing GWAS: $GWAS_ID"

    # Run Swiss Army Knife tool with PICS calculation
      dx run \
      --instance-type=mem3_ssd1_v2_x8 \
      --priority="normal" \
      --name "run_SuSiE_${GWAS_ID}" \
      --tag "SuSiE_calc_gwas" \
      --brief \
      -y \
      swiss-army-knife \
      -iimage_file="project-Gz0KVk8JPx8q9yBzpgXX6y8K:/project-resources/docker_images/tidyverse_docker_w_susie_coloc.tar.gz" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/analysis_scripts/susie_finemapping_analysis/codes/run_susie_for_sig_ld_blocks.R" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/susie_finemapping_analysis/meta_traits/sig_ld_blocks/${GWAS_ID}_meta_out_EUR_significant_ld_blocks.txt" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/meta_EUR_only_results/meta/regenie_format/${GWAS_NAME}" \
      --destination="${OUTPUT_DIR}" \
      -icmd="Rscript run_susie_for_sig_ld_blocks.R \
        --ld_block_file ${GWAS_ID}_meta_out_EUR_significant_ld_blocks.txt \
        --gwas_data_file  ${GWAS_NAME} \
        --gwas_name_prefix ${GWAS_ID}"

    if [ $? -ne 0 ]; then
        echo "Error occurred while processing GWAS: $GWAS_ID"
        exit 1
    fi
    
#done
