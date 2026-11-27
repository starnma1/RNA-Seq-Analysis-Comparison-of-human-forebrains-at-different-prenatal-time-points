#!/usr/bin/env Rscript
# =============================================================================
# Data-pre-processing.R
#
# Prepares the two protein-coding count matrices used downstream:
#   - protein_coding_counts.csv: no low-count filter, for edgeR (spline DE,
#     which filters with filterByExpr) and the CIBERSORTx mixture file
#   - plotting_counts.csv: low-count filtered, for PCA and the heatmap
#
# Steps
#   1. Read the unique-read featureCounts matrix.
#   2. Rename columns to run IDs (ERR...) and order them like sample_metadata.tsv.
#   3. Keep protein-coding genes only (Ensembl gene_biotype), as in
#      Cardoso-Moreira et al. 2019.
#   4. Save that matrix (protein_coding_counts.csv).
#   5. Remove genes with too few reads to be informative
#      (keep genes with >= 10 counts in >= 3 samples) and save
#      plotting_counts.csv.
#
# Inputs
#   quantification/Output/gene_counts_matrix_unique.csv
#   quantification/Output/gene_counts_annotated_unique.csv   (gene names, biotypes)
#   analysis/data/sample_metadata.tsv
#
# Outputs (analysis/data/)
#   protein_coding_counts.csv            protein-coding genes x samples, raw integer
#                                        counts, NO low-count filter (edgeR, CIBERSORTx)
#   plotting_counts.csv                  same, after the low-count filter (PCA, heatmap)
#                                        Both: columns are run IDs (ERR...) in
#                                        sample_metadata.tsv order
#   gene_annotation.tsv                  gene_id, gene_name, gene_biotype, gene_length
#                                        for all protein-coding genes
#   preprocessing_gene_numbers.tsv       how many genes remain after each step
#
# Notes
#   The sex check (04_sex_check.R) deliberately does NOT use this matrix:
#   XIST is a non-coding gene and is removed here.
#   Run from the RStudio project (paths are built with here()).
# =============================================================================

library(here)
library(tidyverse)

data_dir <- here("analysis", "data")
count_dir <- here("quantification", "Output")

# Low-count filter settings
min_count   <- 10   # a gene needs at least this many reads ...
min_samples <- 3    # ... in at least this many samples (~ smallest group size)

# ---- Read inputs -----------------------------------------------------------------
metadata <- read_tsv(file.path(data_dir, "sample_metadata.tsv"))

counts <- as.matrix(read.csv(file.path(count_dir, "gene_counts_matrix_unique.csv"), row.names = 1))

genes <- read_csv(file.path(count_dir, "gene_counts_annotated_unique.csv"),
                  show_col_types = FALSE) %>%
  select(gene_id, gene_name, gene_biotype, gene_length)

# ---- 1. Match columns to the metadata ----------------------------------------
colnames(counts) <- str_extract(colnames(counts), "ERR\\d+")

counts <- counts[, metadata$run]                       # same order as metadata rows

# featureCounts gives whole numbers; store them as integers for DESeq2
storage.mode(counts) <- "integer"

n_all <- nrow(counts)

# ---- 2. Keep protein-coding genes ----------------------------------------------
pc_ids <- genes %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(gene_id)

counts <- counts[rownames(counts) %in% pc_ids, ]
n_protein_coding <- nrow(counts)

# Unfiltered protein-coding matrix (edgeR filters itself with filterByExpr)
write.csv(counts, file.path(data_dir, "protein_coding_counts.csv"))

# ---- 3. Remove low-count genes ---------------------------------------------------
keep   <- rowSums(counts >= min_count) >= min_samples
counts <- counts[keep, ]
n_expressed <- nrow(counts)

# ---- Save ------------------------------------------------------------------------
write.csv(counts, file.path(data_dir, "plotting_counts.csv"))

# Annotation for every protein-coding gene, so names are found for genes in
# either matrix
genes %>%
  filter(gene_biotype == "protein_coding") %>%
  write_tsv(file.path(data_dir, "gene_annotation.tsv"))

tibble(
  step    = c("All genes (featureCounts, GENCODE v50 basic)",
              "Protein-coding",
              paste0("Protein-coding, >= ", min_count, " reads in >= ",
                     min_samples, " samples")),
  n_genes = c(n_all, n_protein_coding, n_expressed)
) %>%
  write_tsv(file.path(data_dir, "preprocessing_gene_numbers.tsv"))
