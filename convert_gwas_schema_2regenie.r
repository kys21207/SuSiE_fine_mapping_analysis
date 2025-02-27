#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(R.utils)

# Function to calculate -log10(p_value)
log10p <- function(p_value) {
  -log10(p_value)
}

# Read input arguments
#args <- commandArgs(trailingOnly = TRUE)
input_file <- "/mnt/project/publically_available_supporting_files/gwas_public_results/autoimmune_step_2_additive_psoriasis_harm_input.h.tsv.gz"
output_file <- "/opt/notebooks/gwas_data/autoimmune_step_2_additive_psoriasis_harm_input.regenie.tsv"

# Read the gzipped input file
data <- fread(cmd = paste("gunzip -c", input_file), sep = "\t", header = TRUE)

# Generate the ID column
data[, ID := paste(chromosome, base_pair_location, other_allele, effect_allele, sep = ":")]

# Calculate LOG10P column
data[, LOG10P := log10p(p_value)]

# Calculate CHISQ column
data[, CHISQ := (beta / standard_error)^2]

# Select and rename columns according to schema1
output_data <- data[, .(CHROM = chromosome, GENPOS = base_pair_location, ID, 
                        ALLELE0 = other_allele, ALLELE1 = effect_allele, 
                        A1FREQ = effect_allele_frequency, N = 426000, BETA = beta, 
                        SE = standard_error, CHISQ, LOG10P, INFO = "ADD")]

# Write the output to a gzipped file
fwrite(output_data, file = gzfile(output_file, "w"), sep = "\t", na = "NA", quote = FALSE)

cat("Conversion complete. Output saved to", output_file, "\n")
