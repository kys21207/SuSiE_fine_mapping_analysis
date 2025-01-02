#!/bin/bash

chmod +x codes/run_pipeline4finemapping_cox.sh 

# Read the file line by line
while IFS= read -r line; do
  # Split the line into arguments
  set -- $line
  # Run the script with the arguments
  ./codes/run_pipeline4finemapping_cox.sh "$1" "$2" "$3" "$4"

echo "Run $2"
  
done < codes/argu_data/argu_data_cox.txt

