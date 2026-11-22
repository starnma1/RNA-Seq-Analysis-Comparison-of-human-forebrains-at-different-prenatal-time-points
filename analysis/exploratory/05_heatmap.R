#!/usr/bin/env Rscript
# =============================================================================
# 05_heatmap.R
#
# Heatmap of the 50 most variable protein-coding genes (TMM-normalised
# log2 CPM, scaled per gene), with samples ordered by age.
#
# Inputs (analysis/data/)
#   processed_count_matrix.csv   genes x samples, raw counts (from Data-pre-processing.R)
#   gene_annotation.tsv          gene_id, gene_name, gene_biotype, gene_length
#   sample_metadata.tsv          run, sample_name, sex, wpc, ... (from 01_samples.R)
#
# Output
#   Results-presentation/Figures/heatmap_top50_variable_genes.png
#   (also shown in the plot pane); prints the gene names in the order the
#   clustering puts them.
#
# Run from the RStudio project (paths are built with here()).
# =============================================================================

library(edgeR)
library(pheatmap)
library(here)

# Settings ----
data_dir   <- here("analysis", "data")
figure_dir <- here("Results-presentation", "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load data ----

counts   <- as.matrix(read.csv(file.path(data_dir, "processed_count_matrix.csv"),
                               row.names = 1, check.names = FALSE))
genes    <- read.delim(file.path(data_dir, "gene_annotation.tsv"))
metadata <- read.delim(file.path(data_dir, "sample_metadata.tsv"))

# Match metadata to count columns, then order samples by age.
sample_index <- match(colnames(counts), metadata$run)
stopifnot(!anyNA(sample_index))
metadata <- metadata[sample_index, ]
age_order <- order(metadata$wpc)
counts <- counts[, age_order]
metadata <- metadata[age_order, ]

# Normalize and select genes by SD before scaling the heatmap.
dge <- DGEList(counts = counts)
dge <- normLibSizes(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)

gene_sd <- apply(log_cpm, 1, sd)
top_genes <- names(sort(gene_sd, decreasing = TRUE))[1:50]
heatmap_data <- log_cpm[top_genes, ]

# Show gene names alongside IDs.
gene_names <- genes$gene_name[match(top_genes, genes$gene_id)]
labels <- gene_names

sample_annotation <- data.frame(
  Weeks = metadata$wpc,
  Sex = metadata$sex,
  row.names = metadata$run
)

hm <- pheatmap(
  heatmap_data,
  main = "Top 50 variable protein-coding genes",
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_method = "ward.D",
  treeheight_row = 0,
  cluster_cols = FALSE,
  annotation_col = sample_annotation,
  labels_row = labels,
  show_colnames = FALSE
)

# Save the heatmap (pheatmap has already drawn it in the plot pane)
png(file.path(figure_dir, "heatmap_top50_variable_genes.png"),
    width = 10, height = 10, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(hm$gtable)
dev.off()

ordered_gene_ids <- rownames(heatmap_data)[hm$tree_row$order]
genes[match(ordered_gene_ids, genes$gene_id), ]$gene_name

