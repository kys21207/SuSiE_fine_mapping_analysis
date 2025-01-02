# prepare_susie_data.R 
# Prepares data for SuSiE for a given LD block.

# Load necessary libraries
library(data.table)

# Get block ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No block ID provided. Usage: Rscript prepare_susie_data.R <BLOCK_ID>")
}
block_id <- args[1]

# Set paths to directories
SUMSTATS_DIR <- "susie_input"
LD_MATRICES_DIR <- "ld_matrices"
OUTPUT_DIR <- "susie_input"

# Load summary statistics
sumstats_file <- file.path(SUMSTATS_DIR, paste0("sumstats_block_", block_id, ".txt"))
sumstats <- fread(sumstats_file)
sumstats <- sumstats[!startsWith(ID, "HGSV")]

#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
colnames(sumstats) <- c("CHROM", "GENPOS", "SNPID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

# Ensure column names are consistent
#colnames(sumstats) <- c("CHROM", "GENPOS", "SNPID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P", "INFO")

# Load LD matrix
ld_matrix_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", block_id, ".ld.gz"))
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
beta <- sumstats$BETA
se <- sumstats$SE
n <- sumstats$N[1]  # Assuming N is the same for all SNPs

# Save prepared data
output_file <- file.path(OUTPUT_DIR, paste0("susie_data_block_", block_id, ".rds"))
saveRDS(list(beta = beta,
             se = se,
             n = n,
             R = ld_matrix,
             snp_ids = snp_ids),
        file = output_file)

