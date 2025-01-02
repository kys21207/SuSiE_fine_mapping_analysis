#!/bin/bash 

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory name and a file names for a GWAS summary statistic."
    echo "and a final result name."
    echo "Usage: $0 GWAS_SUM_dir GWAS_SUM_file_name file_name"
    exit 1
fi

gwasdir=$1
gwasname=$2
finaloutputname=$3
samplesize=$4

chmod +x codes/filter_sumstats_cox.sh 
chmod +x codes/identify_significant_ld_blocks_cox.sh 
chmod +x codes/compute_ld_matrices_for_significant_blocks.sh 
chmod +x codes/prepare_susie_input_for_significant_blocks_cox.sh 
chmod +x codes/run_susie_for_significant_blocks_cox.sh 

./codes/filter_sumstats_cox.sh ${gwasdir} ${gwasname}
./codes/identify_significant_ld_blocks_cox.sh
./codes/compute_ld_matrices_for_significant_blocks.sh 
./codes/prepare_susie_input_for_significant_blocks_cox.sh ${gwasdir} ${gwasname} ${samplesize}
./codes/run_susie_for_significant_blocks_cox.sh ${finaloutputname}

echo "The fine-mapping analysis is completed for ${gwasname} using SuSieR."

# Create the archive for all SuSiE results
tar -czvf "susie_results/${finaloutputname}_susie_results.tar.gz" "susie_results/${finaloutputname}_"*
# Remove all SuSiE results
rm "susie_results/${finaloutputname}_"*.rds "susie_results/${finaloutputname}_"*.txt

