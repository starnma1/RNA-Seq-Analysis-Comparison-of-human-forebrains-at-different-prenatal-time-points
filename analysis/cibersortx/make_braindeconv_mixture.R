#!/usr/bin/env Rscript
# =============================================================================
# make_braindeconv_mixture.R
#
# Mixture file for BrainDeconvShiny (https://voineagulab.shinyapps.io/BrainDeconvShiny/)
# The app wants: .csv or .tsv, Ensembl IDs as rows WITHOUT version numbers and
# unique, sample names as columns, file < 10 MB.
#
# Input   analysis/data/protein_coding_counts.csv   (raw counts, no low-count filter)
#         analysis/data/gene_annotation.tsv         (gene lengths)
# Output  analysis/data/braindeconv_mixture_TPM.csv  TPM, linear scale
# =============================================================================

library(here)
library(tidyverse)

data_dir <- here("analysis", "data")

counts <- as.matrix(read.csv(file.path(data_dir, "protein_coding_counts.csv"),
                             row.names = 1, check.names = FALSE))
genes  <- read_tsv(file.path(data_dir, "gene_annotation.tsv"), show_col_types = FALSE)

len <- genes$gene_length[match(rownames(counts), genes$gene_id)]
ok  <- !is.na(len) & len > 0
counts <- counts[ok, ]; len <- len[ok]

# Counts -> TPM
rpk <- counts / (len / 1000)
tpm <- t(t(rpk) / colSums(rpk)) * 1e6

# Ensembl IDs without version (ENSG...12 -> ENSG...); drop duplicates (e.g. _PAR_Y)
ids  <- sub("\\..*$", "", rownames(tpm))
keep <- !duplicated(ids) & rowSums(tpm) > 0
tpm  <- tpm[keep, ]
rownames(tpm) <- ids[keep]

out <- file.path(data_dir, "braindeconv_mixture_TPM.csv")
write.csv(signif(tpm, 4), out, quote = FALSE)
message("Wrote ", nrow(tpm), " genes x ", ncol(tpm), " samples, ",
        round(file.size(out) / 1e6, 1), " MB")
