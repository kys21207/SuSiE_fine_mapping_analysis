#!/bin/bash 

# Directories and files

EQTL_DATA_PATH="/mnt/project/project-resources/rosmap_brain"
OUTPUT_DIR="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/susie_finemapping_analysis/rosmap_brain"

#for GWAS_FILE in "$GWAS_DATA_PATH"/*.tsv.gz; do
EQTL_FILE="/mnt/project/project-resources/rosmap_brain/celltype-eqtl-sumstats.Ast.tsv.gz"
    EQTL_ID=$(basename "$EQTL_FILE" .tsv.gz)
    EQTL_NAME=$(basename "$EQTL_FILE" )

    # Prepare data in R
    echo "Processing eQTL: $EQTL_ID"

    # Run Swiss Army Knife tool with PICS calculation
      dx run \
      --instance-type=mem3_ssd1_v2_x8 \
      --priority="normal" \
      --name "run_SuSiE_${EQTL_ID}" \
      --tag "SuSiE_calc_eqtl" \
      --brief \
      -y \
      swiss-army-knife \
      -iimage_file="project-Gz0KVk8JPx8q9yBzpgXX6y8K:/project-resources/docker_images/tidyverse_docker_w_susie_coloc.tar.gz" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:analyses_KJ/analysis_scripts/susie_finemapping_analysis/codes/run_susie_for_eqtl_sig_ld_blocks.R" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:project-resources/rosmap_brain/${EQTL_ID}.tsv.gz" \
      -iin="project-Gz0KVk8JPx8q9yBzpgXX6y8K:project-resources/rosmap_brain/rosmap_brain_info.txt" \
      --destination="${OUTPUT_DIR}" \
      -icmd="Rscript run_susie_for_eqtl_sig_ld_blocks.R \
        --eqtl_data_file ${EQTL_ID}.tsv.gz \
        --eqtl_name_prefix  ${EQTL_ID} \
        --eqtl_info_file rosmap_brain_info.txt"

    if [ $? -ne 0 ]; then
        echo "Error occurred while processing eQTL: $EQTL_ID"
        exit 1
    fi
    
#done
