# prepare_susie_data.R 
# Prepares data for SuSiE for a given LD block.

# Load necessary libraries
library(data.table)

# Get block ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) <= 1) {
  stop("No block ID or sample size provided. Usage: Rscript prepare_susie_data.R <BLOCK_ID> <sample szie>")
}
block_id <- args[1]
sample_n <- args[2]

# Set paths to directories
SUMSTATS_DIR <- "susie_input"
LD_MATRICES_DIR <- "ld_matrices"
OUTPUT_DIR <- "susie_input"

# Load summary statistics
sumstats_file <- file.path(SUMSTATS_DIR, paste0("sumstats_block_", block_id, ".txt"))
sumstats <- fread(sumstats_file)
sumstats <- sumstats[!startsWith(ID, "HGSV")]

# Ensure column names are consistent
colnames(sumstats) <- c("CHROM","GENPOS","ALLELE0","ALLELE1","SNPID","MAF","missing.rate","p.value.spa","p.value.norm","Stat","Var","z","A1FREQ")

# Load LD matrix
ld_matrix_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", block_id, ".ld"))
ld_matrix <- as.matrix(fread(ld_matrix_file, header = FALSE))

# Load SNP IDs from BIM file
bim_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", block_id, ".bim"))
bim <- fread(bim_file, header = FALSE)
snp_ids <- bim$V2

# Ensure SNPs are in the same order
sumstats <- sumstats[match(snp_ids, sumstats$SNPID), ]

if (any(is.na(sumstats$SNPID))) {
  cat("Some SNPs in the LD matrix are not in the summary statistics for block", block_id, "\n")
  # Remove SNPs not present in both datasets
  valid_indices <- which(!is.na(sumstats$SNPID))
  sumstats <- sumstats[valid_indices, ]
  ld_matrix <- ld_matrix[valid_indices, valid_indices]
  snp_ids <- snp_ids[valid_indices]
}

# Prepare data for SuSiE
z <- sumstats$z
n <- as.numeric(sample_n)  # Assuming N is the same for all SNPs

# Save prepared data
output_file <- file.path(OUTPUT_DIR, paste0("susie_data_block_", block_id, ".rds"))
saveRDS(list(z = z,
             n = n,
             R = ld_matrix,
             snp_ids = snp_ids),
        file = output_file)

