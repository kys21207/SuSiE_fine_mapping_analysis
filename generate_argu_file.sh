#!/bin/bash

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments supplied. Please provide a directory name for a GWAS summary statistic."
    echo "Usage: $0 GWAS_SUM_dir"
    exit 1
fi

# Define the folder name
folder=$1
name=$2

# Create or clear the argu.txt file
> codes/argu_data/argu_data_${name}.txt

# Loop through each .regenie.tsv.gz file in the folder
#for file in "$folder"/*.regenie.gz; do
for file in "$folder"/*.regenie; do
  # Extract the base name without extension
#  base_name=$(basename "$file" .regenie.gz)
  base_name=$(basename "$file" .regenie)
  
  # Write the folder name, file name, and output file name to argu.txt
  echo "$folder $(basename "$file") $base_name" >> codes/argu_data/argu_data_${name}.txt
done

