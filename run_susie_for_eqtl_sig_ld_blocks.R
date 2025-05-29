# Load necessary libraries
suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(R.utils)
  library(Rfast)
#  library(optparse)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option("--eqtl_data_file", type="character", 
              default=NULL,
              help="Path to eqtl index file", metavar="CHARACTER"),
  make_option("--eqtl_name_prefix", type="character", default=NULL,
              help="gwas file name prefix", metavar="CHARACTER")
  make_option("--eqtl_info_file", type="character", default=NULL,
              help="eqtl sample size file", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#eqtl_data_file="/mnt/project/project-resources/rosmap_brain/celltype-eqtl-sumstats.Ast.tsv.gz"
#eqtl_name_prefix="Ast"
#eqtl_info_file="/mnt/project/project-resources/rosmap_brain/rosmap_brain_info.txt"
#SUMSTATS_FILE   <- eqtl_data_file
#EQTL_ID         <- eqtl_name_prefix

SUMSTATS_FILE_PATH   <- opt$eqtl_data_file
EQTL_ID              <- opt$eqtl_name_prefix
EQTL_INFO_PATH       <-opt$eqtl_info_file
LD_MATRICES_DIR <- "/mnt/project/analyses_KJ/ld_matrices"
LD_BLOCKS_FILE  <- "/mnt/project/analyses_KJ/analysis_scripts/cal_ld_matrices/codes/ld_blocks/ld_blocks_with_ids.bed"

# Load summary statistics
message("Loading summary statistics...")
sumstats <- fread(SUMSTATS_FILE_PATH)

# Load eqtl sample size
message("Loading eqtl sample size...")
eqtl_info <- read.table(EQTL_INFO_PATH, header=T, sep="")
N <- eqtl_info$sample_size[eqtl_info$study_id == EQTL_ID]

# Load LD blocks
message("Loading LD blocks...")
ld_blocks <- fread(LD_BLOCKS_FILE, header = FALSE)
colnames(ld_blocks) <- c("chr", "start", "end", "block_id")

# Add ID column to SNPs
message("Creating SNP IDs...")
sumstats[, ID := paste0(gsub("chr","",chr38),":",pos38,":",REF,":",ALT)]

# OPTIMIZED SECTION: Assign SNPs to LD blocks using foverlaps
message("Preparing data for efficient block assignment...")

# Ensure both tables are data.tables
setDT(sumstats)
setDT(ld_blocks)

# Make sure chromosome columns are the same type
ld_blocks[, chr := as.character(chr)]
sumstats[, chr38 := as.character(chr38)]

# Prepare data for overlap join
# For SNPs, start and end are the same (the position)
message("Setting up interval join...")
snps_intervals <- sumstats[, .(
  chr = chr38,
  start = pos38,
  end = pos38,
  row_id = .I
)]

# Rename ld_blocks columns to match what foverlaps expects
setnames(ld_blocks, c("start", "end"), c("start_block", "end_block"))

# Set keys for the overlap join
setkey(ld_blocks, chr, start_block, end_block)
setkey(snps_intervals, chr, start, end)

# Perform the overlap join
message("Performing interval join...")
block_assignments <- foverlaps(snps_intervals, ld_blocks, type="within")

# Update the original sumstats table with block_id
message("Updating summary statistics with block IDs...")
sumstats[, block_id := NA_integer_]
sumstats[block_assignments[!is.na(block_id), row_id], block_id := block_assignments[!is.na(block_id), block_id]]

# Identify the primary block for each gene (block with the most significant SNP)
message("Determining primary LD block for each gene...")
gene_block_map <- sumstats[!is.na(block_id), .(
  min_pvalue = min(pvalue, na.rm = TRUE),
  block_id = block_id[which.min(pvalue)]
), by = gene_id]

# Filter genes with significant SNPs (p < 1e-5)
gene_block_map <- gene_block_map[min_pvalue <= 1e-5]

# Group genes by block_id for efficient processing
gene_groups <- split(gene_block_map$gene_id, gene_block_map$block_id)

message(sprintf("Processing %d LD blocks containing significant genes...", length(gene_groups)))

# Create a combined data table to accumulate all credible sets
all_credible_sets <- data.table()

# Process each block
for (grp_block_id in names(gene_groups)) {
  block_genes <- gene_groups[[grp_block_id]]
  message(sprintf("Processing block %s with %d genes", grp_block_id, length(block_genes)))
  
  # Skip small gene groups
  if (length(block_genes) == 0) {
    message(sprintf("No genes to process in block %s. Skipping.", grp_block_id))
    next
  }
  
  # Load LD matrix just once per block
  ld_matrix_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", grp_block_id, ".ld.gz"))
  vars_file <- file.path(LD_MATRICES_DIR, paste0("ld_matrix_block_", grp_block_id, ".vars"))
  
  if (!file.exists(ld_matrix_file) || !file.exists(vars_file)) {
    message(sprintf("LD matrix files for block %s not found. Skipping.", grp_block_id))
    next
  }
  
  # Load LD matrix and SNP IDs
  message("Loading LD matrix...")
  ld_matrix <- as.matrix(fread(ld_matrix_file, header = FALSE))
  bim <- fread(vars_file, header = FALSE)
  colnames(bim) <- c("ID") # Adjust based on actual columns
  snp_ids <- bim$ID
  
  # Clean LD matrix
  ld_matrix[is.nan(ld_matrix)] <- 0
  diag(ld_matrix) <- 1  # Ensure diagonal is 1
  
  # Get block-specific SNPs
  block_snps <- sumstats[block_id == as.integer(grp_block_id)]
  
  # Process each gene in the block
  for (blk_gene_id in block_genes) {
    message(sprintf("Processing gene %s in block %s", blk_gene_id, grp_block_id))
    
    # Filter SNPs for this gene
    gene_sumstats <- sumstats[gene_id == blk_gene_id & block_id == as.integer(grp_block_id)]
    
    # Output file paths
    susie_output_file <- file.path(paste0(EQTL_ID, "_susie_results_gene_", blk_gene_id, ".rds"))
#    credible_sets_file <- file.path(paste0(EQTL_ID, "_credible_sets_gene_", blk_gene_id, ".txt"))
    
    # Skip if results already exist
    if (file.exists(susie_output_file) && file.exists(credible_sets_file)) {
      message(sprintf("Results for gene %s already exist. Skipping.", blk_gene_id))
      next
    }
    
    # Check if we have enough SNPs
    if (nrow(gene_sumstats) < 10) {
      message(sprintf("Gene %s has only %d SNPs. Skipping.", blk_gene_id, nrow(gene_sumstats)))
      next
    }
    
    tryCatch({
      # Match SNPs between LD matrix and gene summary statistics
      message("Matching SNPs between LD matrix and summary statistics...")
      matched_indices <- match(snp_ids, gene_sumstats$ID)
      gene_sumstats_matched <- gene_sumstats[matched_indices[!is.na(matched_indices)]]
      
      # Get corresponding indices in the LD matrix
      valid_indices <- which(!is.na(matched_indices))
      if (length(valid_indices) < 10) {
        message(sprintf("After matching, gene %s has only %d valid SNPs. Skipping.", 
                       blk_gene_id, length(valid_indices)))
        next
      }
      
      # Subset LD matrix for this gene's SNPs
      gene_ld_matrix <- ld_matrix[valid_indices, valid_indices]
      gene_snps <- snp_ids[valid_indices]
      
      # Prepare data for SuSiE
      message("Running SuSiE...")
      beta <- gene_sumstats_matched$beta
      se <- gene_sumstats_matched$se
      n <- gene_sumstats_matched$N  # Assuming N is the same for all SNPs
      
      # Run SuSiE
      susie_fit <- susie_rss(bhat = beta, shat = se, R = gene_ld_matrix, n = n,
                           L = 10,  # Maximum number of causal signals
                           estimate_residual_variance = FALSE,
                           estimate_prior_variance = TRUE)
      
      # Save SuSiE results
      saveRDS(susie_fit, susie_output_file)
      message(sprintf("SuSiE results saved to %s", susie_output_file))
      
      # Extract credible sets
      if (!is.null(susie_fit$sets$cs) && length(susie_fit$sets$cs) > 0) {
        pips <- susie_fit$pip
        
        credible_sets <- data.table()
        cs_indices <- susie_fit$sets$cs
        
        for (j in seq_along(cs_indices)) {
          cs_snps <- gene_snps[cs_indices[[j]]]
          cs_pip <- pips[cs_indices[[j]]]
          
          cs_dt <- data.table(
            Block = grp_block_id,
            Gene = blk_gene_id,
            CS_Number = j,
            SNP = cs_snps,
            PIP = cs_pip
          )
          
          credible_sets <- rbind(credible_sets, cs_dt)
        }
      
        # Append this block's credible sets to the combined results
        all_credible_sets <- rbind(all_credible_sets, credible_sets)
      } else {
        message(sprintf("No credible sets found for gene %s", blk_gene_id))
      }
      
    }, error = function(e) {
      message(sprintf("Error processing gene %s: %s", blk_gene_id, e$message))
    })
  }
}

# Save all credible sets outside the loop
# Create a single combined output file with all results
combined_output_file <- paste0(EQTL_ID, "_all_credible_sets.txt")
fwrite(all_credible_sets, combined_output_file, sep = "\t")
message(sprintf("All credible sets saved to %s", combined_output_file))

message("All blocks processed.")
