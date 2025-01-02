# process_ld_matrix.R 
# Reads the LD matrix file, converts it to single precision,
# saves it as a compressed RDS file.

args <- commandArgs(trailingOnly = TRUE)
ld_file <- args[1]          # Path to the LD matrix file (from PLINK)
output_file <- args[2]      # Path to the output RDS file (compressed)

# Load the LD matrix
ld_matrix <- as.matrix(read.table(ld_file, header = FALSE))

# Convert to single-precision float
ld_matrix_single <- matrix(as.single(ld_matrix), nrow = nrow(ld_matrix))

# Save as compressed RDS file
saveRDS(ld_matrix_single, gzfile(output_file))


