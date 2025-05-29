# Run data for SuSiE for given LD blocks.

# Load necessary libraries
suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(R.utils)
  library(Rfast)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option("--ld_block_file", type="character", default=NULL, 
              help="Path to the ld block info file", metavar="CHARACTER"),
  make_option("--gwas_data_file", type="character", 
              default=NULL,
              help="Path to eqtl index file", metavar="CHARACTER"),
  make_option("--gwas_name_prefix", type="character", default=NULL,
              help="gwas file name prefix", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

SUMSTATS_FILE   <- opt$gwas_data_file
BLOCK_FILE      <- opt$ld_block_file
GWAS_ID         <- opt$gwas_name_prefix

#BLOCK_FILE <- "/opt/notebooks/tmp/G6_MS_meta_out_EUR_significant_ld_blocks.txt"
#SUMSTATS_FILE <- "/mnt/project/analyses_KJ/meta_EUR_only_results/meta/regenie_format/G6_MS_meta_out_EUR.tsv.gz"
#GWAS_ID <-"G6_MS"
LD_MATRICES_DIR <- "/mnt/project/analyses_KJ/ld_matrices"

# Load block information
message("Loading block information...")
block_data <- fread(BLOCK_FILE, header = FALSE, col.names = c("CHR", "START", "END", "BLOCK_ID"))

# Load summary statistics once
message("Loading summary statistics...")
sumstats <- fread(SUMSTATS_FILE)

# Make sure required columns exist
required_cols <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE")
if (!all(required_cols %in% colnames(sumstats))) {
  stop("Summary statistics must contain the following columns: ", paste(required_cols, collapse = ", "))
}

# Create a combined data table to accumulate all credible sets
all_credible_sets <- data.table()

# Process each block
for (i in 1:nrow(block_data)) {
  block_row <- block_data[i,]
  BLOCK_ID <- block_row$BLOCK_ID
  chr <- as.numeric(gsub("chr","",block_row$CHR))
  start_pos <- block_row$START
  end_pos <- block_row$END
  
  message(sprintf("Processing block %s (chr%s:%s-%s)...", BLOCK_ID, chr, start_pos, end_pos))
  
  # Output file paths
  susie_output_file <- file.path(paste0(GWAS_ID, "_susie_results_block_", BLOCK_ID, ".rds"))
  credible_sets_file <- file.path(paste0(GWAS_ID, "_credible_sets_block_", BLOCK_ID, ".txt"))
  
  # Skip if results already exist
  if (file.exists(susie_output_file) && file.exists(credible_sets_file)) {
    message(sprintf("Results for block %s already exist. Skipping.", BLOCK_ID))
    next
  }
  
  # Filter summary statistics for current block
  block_sumstats <- sumstats[CHROM == chr & 
                             GENPOS >= start_pos & 
                             GENPOS <= end_pos]
  
  # Remove non-standard variants if SNPID starts with "HGSV"
  block_sumstats <- block_sumstats[!startsWith(ID, "HGSV")]
  
  # Check if we have enough SNPs
  if (nrow(block_sumstats) < 10) {
    message(sprintf("Block %s has only %d SNPs. Skipping.", BLOCK_ID, nrow(block_sumstats)))
    next
  }
  
  # Load LD matrix
  ld_matrix_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", BLOCK_ID, ".ld.gz"))
  vars_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", BLOCK_ID, ".vars"))
  
  if (!file.exists(ld_matrix_file) || !file.exists(vars_file)) {
    message(sprintf("LD matrix files for block %s not found. Skipping.", BLOCK_ID))
    next
  }
  
  tryCatch({
    # Load LD matrix and SNP IDs
    message("Loading LD matrix...")
    ld_matrix <- as.matrix(fread(ld_matrix_file, header = FALSE))
    bim <- fread(vars_file, header = FALSE)
    snp_ids <- bim$V1
    
    # Clean LD matrix
    ld_matrix[is.nan(ld_matrix)] <- 0
    diag(ld_matrix) <- 1  # Ensure diagonal is 1
    
    # Match SNPs between LD matrix and summary statistics
    message("Matching SNPs between LD matrix and summary statistics...")
    block_sumstats <- block_sumstats[match(snp_ids, ID)]
    
    # Remove SNPs not present in both datasets
    valid_indices <- which(!is.na(block_sumstats$ID))
    if (length(valid_indices) < 10) {
      message(sprintf("After matching, block %s has only %d valid SNPs. Skipping.", 
                     BLOCK_ID, length(valid_indices)))
      next
    }
    
    block_sumstats <- block_sumstats[valid_indices]
    ld_matrix <- ld_matrix[valid_indices, valid_indices]
    snp_ids <- snp_ids[valid_indices]
    
    # Prepare data for SuSiE
    message("Running SuSiE...")
    beta <- block_sumstats$BETA
    se <- block_sumstats$SE
    n <- block_sumstats$N[1]  # Assuming N is the same for all SNPs
    
    # Run SuSiE with proper error handling
    susie_fit <- susie_rss(bhat = beta, shat = se, R = ld_matrix, n = n,
                          L = 10,  # Maximum number of causal signals
                          estimate_residual_variance = FALSE,
                          estimate_prior_variance = TRUE)
    
    # Save SuSiE results
    saveRDS(susie_fit, susie_output_file)
    message(sprintf("SuSiE results saved to %s", susie_output_file))
    
    # Extract credible sets and accumulate them
    if (!is.null(susie_fit$sets$cs) && length(susie_fit$sets$cs) > 0) {
      pips <- susie_fit$pip
      
      credible_sets <- data.table()
      cs_indices <- susie_fit$sets$cs
      
      for (j in seq_along(cs_indices)) {
        cs_snps <- snp_ids[cs_indices[[j]]]
        cs_pip <- pips[cs_indices[[j]]]
        
        cs_dt <- data.table(
          Block = BLOCK_ID,
          CS_Number = j,
          SNP = cs_snps,
          PIP = cs_pip,
          OutputFile = credible_sets_file  # Store the target output file
        )
        
        credible_sets <- rbind(credible_sets, cs_dt)
      }
      
      # Append this block's credible sets to the combined results
      all_credible_sets <- rbind(all_credible_sets, credible_sets)
    } else {
      message(sprintf("No credible sets found for block %s", BLOCK_ID))
    }
    
  }, error = function(e) {
    message(sprintf("Error processing block %s: %s", BLOCK_ID, e$message))
  })
}

# Save all credible sets outside the loop
# Create a single combined output file with all results
combined_output_file <- paste0(GWAS_ID, "_all_credible_sets.txt")
fwrite(all_credible_sets[, .(Block, CS_Number, SNP, PIP, OutputFile)], combined_output_file, sep = "\t")
message(sprintf("All credible sets saved to %s", combined_output_file))

# Create a script that will split the combined file into individual files
#split_script_file <- paste0(GWAS_ID, "_split_credible_sets.R")
#split_script <- c(
#  "# Auto-generated script to split credible sets into individual files",
#  "library(data.table)",
#  paste0("combined_data <- fread('", combined_output_file, "')"),
#  "output_files <- unique(combined_data$OutputFile)",
#  "for (file in output_files) {",
#  "  file_data <- combined_data[OutputFile == file, .(Block, CS_Number, SNP, PIP)]",
#  "  fwrite(file_data, file, sep = '\\t')",
#  "  cat(sprintf('Wrote %d rows to %s\\n', nrow(file_data), file))",
#  "}"
#)

# Write the split script but don't execute it
#writeLines(split_script, split_script_file)
#message(sprintf("Created split script at %s", split_script_file))
#message("To split the combined file into individual files, run the generated R script.")

message("All blocks processed.")
