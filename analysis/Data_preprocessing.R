#!/usr/bin/env Rscript
# =============================================================================
# Data-pre-processing.R
#
# Prepares the gene count matrix for all downstream analyses (exploratory,
# DESeq2, GSEA).
#
# Steps
#   1. Read the unique-read featureCounts matrix.
#   2. Rename columns to run IDs (ERR...) and order them like sample_metadata.tsv.
#   3. Keep protein-coding genes only (Ensembl gene_biotype), as in
#      Cardoso-Moreira et al. 2019.
#   4. Remove genes with too few reads to be informative
#      (keep genes with >= 10 counts in >= 3 samples).
#
# Inputs
#   quantification/Output/gene_counts_matrix_unique.csv
#   quantification/Output/gene_counts_annotated_unique.csv   (gene names, biotypes)
#   analysis/data/sample_metadata.tsv
#
# Outputs (analysis/data/)
#   counts_protein_coding_filtered.csv   genes x samples, raw integer counts (DESeq2 input)
#   gene_annotation.tsv                  gene_id, gene_name, gene_biotype, gene_length
#                                        for the genes kept
#   preprocessing_gene_numbers.tsv       how many genes remain after each step
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

# ---- 2. Keep protein-coding genes ----------------------------------------------
pc_ids <- genes %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(gene_id)

counts <- counts[rownames(counts) %in% pc_ids, ]


# ---- 3. Remove low-count genes ---------------------------------------------------
keep   <- rowSums(counts >= min_count) >= min_samples
counts <- counts[keep, ]

# ---- Save ------------------------------------------------------------------------
write.csv(counts, file.path(data_dir, "processed_count_matrix.csv"))

agenes %>%
  filter(gene_id %in% rownames(counts)) %>%
  write_tsv(file.path(data_dir, "gene_annotation.tsv"))

