library(edgeR)
library(pheatmap)
library(here)

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

ordered_gene_ids <- rownames(heatmap_data)[hm$tree_row$order]
genes[match(ordered_gene_ids, genes$gene_id), ]$gene_name

