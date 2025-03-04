# SuSiE_fine_mapping_analysis
## Workflow
                   ┌────────────────────────────────────────────────┐ 
                   │ 1. GWAS Summary Stats + LD Blocks Definition   │
                   │    - Identify all LD blocks (e.g., from        │
                   │      Berisa & Pickrell or other sources)       │
                   └─────────┬──────────────────────────────────────┘
                             │
                             ▼
                   ┌────────────────────────────────────────────────┐
                   │ 2. For Each LD Block:                          │
                   │    a) Subset Summary Stats to SNPs in Block    │
                   │    b) Subset Reference Genotype to SNPs in     │
                   │       Block                                    │
                   │    c) Compute Local LD Matrix                  │
                   └─────────┬──────────────────────────────────────┘
                             │
                             ▼
                   ┌────────────────────────────────────────────────┐
                   │ 3. Prepare SuSiE Inputs for Each Block:        │
                   │    - bhat (Effect Sizes)                       │
                   │    - shat (Standard Errors)                    │
                   │    - R (LD Matrix for This Block)              │
                   │    - n (GWAS Sample Size)                      │
                   └─────────┬──────────────────────────────────────┘
                             │
                             ▼
                   ┌────────────────────────────────────────────────┐
                   │ 4. Run SuSiE per Block:                        │
                   │    susie_rss(bhat, shat, R, n, ...)            │
                   │    - Obtain Posterior Inclusion Probabilities  │
                   │      (PIPs)                                    │
                   │    - Extract Credible Sets (CS)                │
                   └─────────┬──────────────────────────────────────┘
                             │
                             ▼
                   ┌────────────────────────────────────────────────┐
                   │ 5. Combine Results Across Blocks:              │
                   │    - Merge CS from Each Block                  │
                   │    - Summarize Final PIPs & CS                 │
                   └─────────┬──────────────────────────────────────┘
                             │
                             ▼
                   ┌────────────────────────────────────────────────┐
                   │ 6. Post-Processing & Interpretation            │
                   │    - Annotate CS Variants (e.g., VEP/ANNOVAR)  │
                   │    - v2g Mapping (Colocalization, Chromatin)   │
                   │    - Visualization (PIPs, LocusZoom, etc.)     │ 
                   └────────────────────────────────────────────────┘

## Detailed Explanation of Each Step
1.	GWAS Summary Stats + LD Blocks Definition <br>
o	GWAS Summary Statistics: Contains effect sizes (betas), standard errors, p‐values, and SNP identifiers.  <br>
o	LD Blocks: Defined sets of contiguous SNPs expected to share high linkage disequilibrium.  <br>
o	Tools to Define Blocks: Berisa & Pickrell’s blocks, EUR LD blocks from ijmacdon/LDblocks_GRCh38, or any published block definitions.  <br>

2.	For Each LD Block  <br>
o	Subset GWAS Summary Stats: Extract only the SNPs residing in that block.  <br>
o	Subset Reference Genotype: Select the same set of SNPs from reference data (e.g., 1000 Genomes, UK Biobank).  <br>
o	Compute Local LD Matrix: Typically done with PLINK (--r / --r2) or other LD‐computing routines for the SNP set in that block.  <br>

3.	Prepare SuSiE Inputs for Each Block  <br>
o	bhat: Vector of effect sizes (betas) from your GWAS summary for the block’s SNPs.  <br>
o	shat: Vector of standard errors for those betas.  <br>
o	R: LD matrix computed in Step 2.  <br>
o	n: GWAS sample size.  <br>
o	Note: Align alleles carefully (check reference / alternate alleles) and ensure consistent SNP ordering between bhat/shat and the rows/columns in R.  <br>

4.	Run SuSiE per Block  <br>
o	Function: susie_rss(bhat, shat, R, n, ...).  <br>
o	Outputs:  <br>
1.	Posterior Inclusion Probabilities (PIPs) for each SNP in the block.  <br>
2.	Credible Sets (CS): Groups of SNPs that collectively have a high probability (e.g., 95%) of containing the causal variant(s).  <br>

5.	Combine Results Across Blocks  <br>
o	Merge Credible Sets: Gather the CS from each block into a single summary.  <br>
o	Summarize PIPs: Focus on SNPs with high PIPs across blocks or examine block‐by‐block.  <br>
6.	Post‐Processing & Interpretation  <br>
o	Variant Annotation: Use VEP, ANNOVAR, or other annotation tools.  <br>
o	v2g Mapping: Link fine‐mapped variants to target genes via eQTL colocalization, chromatin interaction data, etc.  <br>
o	Visualization: Plot PIPs, create LocusZoom plots, or highlight variants in regional association plots.  <br>
o	Biological Follow‐Up: Investigate identified credible set variants in labs or cross‐reference literature. <br>

