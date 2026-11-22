#!/usr/bin/env Rscript
# =============================================================================
# DE_spline_enrichment.R
#
# Over-representation analysis (ORA) of trajectory clusters from
# DE_spline_interpretation.R using fgsea::fora() and MSigDB pathways.
#
# Inputs
#   Results-presentation/analysis-output/spline_DE_clusters.csv
#   Results-presentation/analysis-output/spline_fit.rds
#
# Outputs
#   Results-presentation/analysis-output/fora_hallmark_results.csv
#   Results-presentation/analysis-output/fora_gobp_results.csv
#   Results-presentation/Figures/fora_hallmark_dotplot.png
#   Results-presentation/Figures/fora_gobp_dotplot.png
#   Results-presentation/Figures/fora_gobp_bars_by_cluster.png
#   Results-presentation/Figures/fora_top_pathway_trajectories.png
# =============================================================================

library(tidyverse)
library(data.table)
library(fgsea)
library(msigdbr)
library(splines)
library(here)


# Settings ----
padj_cutoff <- 0.05
min_size    <- 15
max_size    <- 500

result_dir <- here("Results-presentation", "analysis-output")
figure_dir <- here("Results-presentation", "Figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# 1. Load data and prepare gene identifiers ----
cluster_tab <- read_csv(
  file.path(result_dir, "spline_DE_clusters.csv"),
  show_col_types = FALSE
)
sp <- readRDS(file.path(result_dir, "spline_fit.rds"))

# Remove Ensembl version suffixes (e.g. ENSG00000141510.1 -> ENSG00000141510)
cluster_tab <- cluster_tab %>%
  mutate(ensembl_clean = sub("\\..*", "", gene_id))

universe_genes <- unique(sub("\\..*", "", rownames(sp$y)))
gene_clusters  <- lapply(
  split(cluster_tab$ensembl_clean, cluster_tab$cluster),
  unique
)

all_cluster_ids <- as.character(sort(unique(cluster_tab$cluster)))
message("Loaded ", length(universe_genes), " background genes across ",
        length(all_cluster_ids), " clusters.")


# 2. Fetch MSigDB pathways and build ID-to-Symbol lookup ----
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(!is.na(ensembl_gene), ensembl_gene != "")

gobp_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(!is.na(ensembl_gene), ensembl_gene != "")

hallmark_pathways <- lapply(split(hallmark_df$ensembl_gene, hallmark_df$gs_name), unique)
gobp_pathways     <- lapply(split(gobp_df$ensembl_gene, gobp_df$gs_name), unique)

# Map unversioned Ensembl IDs to gene symbols for readable tables and plots
symbol_lookup <- bind_rows(
  hallmark_df %>% select(ensembl_gene, gene_symbol),
  gobp_df %>% select(ensembl_gene, gene_symbol)
) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(ensembl_gene, .keep_all = TRUE)

id_to_symbol <- setNames(symbol_lookup$gene_symbol, symbol_lookup$ensembl_gene)


# 3. Helper function to run fora() across all clusters ----
run_fora_all <- function(pathways, clusters, universe, id_map, min_sz, max_sz) {
  res_list <- lapply(names(clusters), function(cl) {
    res <- fora(
      pathways = pathways,
      genes    = clusters[[cl]],
      universe = universe,
      minSize  = min_sz,
      maxSize  = max_sz
    )
    res[, cluster := cl]
    return(res)
  })
  
  rbindlist(res_list) %>%
    as_tibble() %>%
    mutate(
      cluster    = factor(cluster, levels = names(clusters)),
      gene_ratio = overlap / size,
      overlap_symbols = vapply(overlapGenes, function(ids) {
        syms <- ifelse(ids %in% names(id_map), id_map[ids], ids)
        paste(syms, collapse = ";")
      }, character(1))
    ) %>%
    arrange(cluster, padj, pval)
}

hallmark_res <- run_fora_all(
  hallmark_pathways, gene_clusters, universe_genes, id_to_symbol, min_size, max_size
)
gobp_res <- run_fora_all(
  gobp_pathways, gene_clusters, universe_genes, id_to_symbol, min_size, max_size
)

# Save annotated CSV tables (drop raw list-column before writing)
write_csv(
  hallmark_res %>% select(-overlapGenes),
  file.path(result_dir, "fora_hallmark_results.csv")
)
write_csv(
  gobp_res %>% select(-overlapGenes),
  file.path(result_dir, "fora_gobp_results.csv")
)


# 4. Plot 1: Hallmark dot plot across all clusters ----
sig_hallmark <- hallmark_res %>%
  filter(padj < padj_cutoff) %>%
  group_by(cluster) %>%
  slice_min(order_by = padj, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    pathway_label = pathway %>%
      sub("^HALLMARK_", "", .) %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 36)
  )

p_hallmark_dot <- ggplot(
  sig_hallmark,
  aes(x = cluster, y = reorder(pathway_label, -padj))
) +
  geom_point(aes(size = gene_ratio, color = -log10(padj))) +
  scale_x_discrete(drop = FALSE) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title    = "Enriched Hallmark pathways by trajectory cluster",
    subtitle = paste0("Top 5 terms per cluster (FDR < ", padj_cutoff, ")"),
    x        = "Trajectory cluster (early -> late peak)",
    y        = NULL,
    size     = "Gene ratio",
    color    = "-log10(FDR)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

print(p_hallmark_dot)
ggsave(
  file.path(figure_dir, "fora_hallmark_dotplot.png"),
  p_hallmark_dot, width = 8.5, height = 6, dpi = 300, bg = "white"
)


# 5. Plot 2: GO:BP dot plot across all clusters ----
# Select the top 5 pathways per cluster by raw p-value
sig_gobp_top <- gobp_res %>%
  group_by(cluster) %>%
  slice_min(order_by = pval, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    pathway_label = pathway %>%
      sub("^GOBP_", "", .) %>%
      str_replace_all("_", " ") %>%
      str_to_sentence() %>%
      str_wrap(width = 42)
  )

p_gobp_dot <- ggplot(
  sig_gobp_top,
  aes(x = cluster, y = reorder(pathway_label, -padj))
) +
  geom_point(aes(size = gene_ratio, color = -log10(padj))) +
  scale_x_discrete(drop = FALSE) +
  scale_color_viridis_c(option = "viridis") +
  labs(
    title    = "Enriched GO Biological Processes by trajectory cluster",
    subtitle = paste0("Top 5 terms per cluster (FDR < ", padj_cutoff, ")"),
    x        = "Trajectory cluster (early -> late peak)",
    y        = NULL,
    size     = "Gene ratio",
    color    = "-log10(FDR)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

print(p_gobp_dot)
ggsave(
  file.path(figure_dir, "fora_gobp_dotplot.png"),
  p_gobp_dot, width = 9.5, height = 8.5, dpi = 300, bg = "white"
)


# 6. Plot 3: Faceted fold-enrichment bar chart per cluster ----
p_gobp_bars <- ggplot(
  sig_gobp_top,
  aes(x = foldEnrichment, y = reorder(pathway_label, foldEnrichment), fill = -log10(padj))
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = paste0(overlap, "/", size)),
    hjust = -0.15, size = 3, color = "grey20"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  facet_wrap(~ paste("Cluster", cluster), scales = "free_y", ncol = 2) +
  labs(
    title    = "Top enriched GO:BP terms inside each cluster",
    subtitle = "Bar labels show overlap genes / pathway size",
    x        = "Fold enrichment",
    y        = NULL,
    fill     = "-log10(FDR)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

print(p_gobp_bars)
ggsave(
  file.path(figure_dir, "fora_gobp_bars_by_cluster.png"),
  p_gobp_bars, width = 12, height = 9, dpi = 300, bg = "white"
)


# 7. Plot 4: Fitted gene curves for the #1 GO:BP term in each cluster ----
# Rebuild scaled trajectories on the 100-point wpc grid from spline_fit.rds
wpc_grid    <- seq(min(sp$meta$wpc), max(sp$meta$wpc), length.out = 100)
grid_basis  <- predict(sp$age_basis, newx = wpc_grid)
grid_design <- matrix(
  0, nrow = length(wpc_grid), ncol = ncol(sp$design),
  dimnames = list(NULL, colnames(sp$design))
)
grid_design[, "Intercept"]       <- 1
grid_design[, sp$spline_coefs]   <- grid_basis

grid_logcpm <- (sp$fit$coefficients %*% t(grid_design) + log(1e6)) / log(2)
rownames(grid_logcpm) <- sub("\\..*", "", rownames(grid_logcpm))
grid_logcpm <- grid_logcpm[!duplicated(rownames(grid_logcpm)), , drop = FALSE]
colnames(grid_logcpm) <- wpc_grid

z_mat <- t(scale(t(grid_logcpm)))

# Extract the single most significant GO:BP term per cluster
top1_per_cluster <- gobp_res %>%
  filter(padj < padj_cutoff) %>%
  group_by(cluster) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  ungroup()

traj_pathway_df <- pmap_dfr(top1_per_cluster, function(cluster, pathway, overlapGenes, ...) {
  valid_ids <- intersect(overlapGenes, rownames(z_mat))
  clean_title <- pathway %>%
    sub("^GOBP_", "", .) %>%
    str_replace_all("_", " ") %>%
    str_to_sentence() %>%
    str_wrap(width = 34)
  
  as_tibble(z_mat[valid_ids, , drop = FALSE], rownames = "ensembl_clean") %>%
    pivot_longer(-ensembl_clean, names_to = "wpc", values_to = "z") %>%
    mutate(
      wpc         = as.numeric(wpc),
      panel_title = paste0("Cluster ", cluster, ":\n", clean_title)
    )
})

mean_pathway_df <- traj_pathway_df %>%
  group_by(panel_title, wpc) %>%
  summarise(z = mean(z), .groups = "drop")

p_pathway_traj <- ggplot(traj_pathway_df, aes(x = wpc, y = z)) +
  geom_line(aes(group = ensembl_clean), color = "steelblue", alpha = 0.15, linewidth = 0.4) +
  geom_line(data = mean_pathway_df, color = "#B2182B", linewidth = 1.2) +
  facet_wrap(~ panel_title, ncol = 3) +
  labs(
    title    = "Expression trajectories of genes in each cluster's top GO:BP term",
    subtitle = "Blue lines = individual pathway genes; red line = pathway mean",
    x        = "Age (wpc)",
    y        = "Scaled expression (z-score)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

print(p_pathway_traj)
ggsave(
  file.path(figure_dir, "fora_top_pathway_trajectories.png"),
  p_pathway_traj, width = 11, height = 7.5, dpi = 300, bg = "white"
)