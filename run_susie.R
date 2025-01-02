# run_susie.R 
# Runs SuSiE for a given LD block.

# Load necessary libraries
library(susieR)
library(Rfast)
library(data.table)

# Get block ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No block ID provided. Usage: Rscript run_susie.R <BLOCK_ID>")
}
block_id <- args[1]
output_name <- args[2]

# Set paths to directories
INPUT_DIR <- "susie_input"
OUTPUT_DIR <- "susie_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Load prepared data
input_file <- file.path(INPUT_DIR, paste0("susie_data_block_", block_id, ".rds"))
if (!file.exists(input_file)) {
  stop(paste("Input data file for block", block_id, "not found at", input_file))
}
data <- readRDS(input_file)

# Extract data components
beta <- data$beta
se <- data$se
R <- data$R
n <- data$n
snp_ids <- data$snp_ids

# Symmetrize the R matrix
# R <- (R + t(R)) / 2

# Run SuSiE
susie_fit <- susie_rss(bhat = beta, shat = se, R = R, n = n)

# Save SuSiE results
output_file <- file.path(OUTPUT_DIR, paste0(output_name,"_susie_results_block_", block_id, ".rds"))
saveRDS(susie_fit, output_file)

# Optionally, extract credible sets and save them
if (!is.null(susie_fit$sets$cs)) {
  pips <- susie_fit$pip

  credible_sets <- data.frame()
  cs_indices <- susie_fit$sets$cs

  for (i in seq_along(cs_indices)) {
    cs_snps <- snp_ids[cs_indices[[i]]]
    cs_pip <- pips[cs_indices[[i]]]

    cs_df <- data.frame(
      Block = block_id,
      CS_Number = i,
      SNP = cs_snps,
      PIP = cs_pip,
      stringsAsFactors = FALSE
    )

    credible_sets <- rbind(credible_sets, cs_df)
  }

  # Save credible sets
  credible_sets_file <- file.path(OUTPUT_DIR, paste0(output_name,"_credible_sets_block_", block_id, ".txt"))
  write.table(credible_sets, credible_sets_file, sep = "\t", row.names = FALSE, quote = FALSE)
}





