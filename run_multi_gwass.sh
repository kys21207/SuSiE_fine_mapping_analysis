#!/bin/bash

# Argument parsing
while getopts "d:n:" flag; do
    case "${flag}" in
        d) GWAS_SUM_DIR=${OPTARG};;
        n) DATA_GRP_NAME=${OPTARG};;
    esac
done

if [[ -z "${GWAS_SUM_DIR}" ]]; then
    echo "GWAS result directory argument is required."
    exit 1
fi

if [[ -z "${DATA_GRP_NAME}" ]]; then
    echo "Data group argument is required."
    exit 1
fi

chmod +x codes/generate_argu_file.sh
chmod +x codes/run_pipeline4finemapping.sh

# Generate argu_data file from a result folder
#./codes/generate_argu_file.sh $GWAS_SUM_DIR $DATA_GRP_NAME

# Read the file line by line
while IFS= read -r line; do
  # Split the line into arguments
  set -- $line
  # Run the script with the arguments
  ./codes/run_pipeline4finemapping.sh "$1" "$2" "$3"

echo "Run $2"
  
done < codes/argu_data/argu_data_${DATA_GRP_NAME}.txt

