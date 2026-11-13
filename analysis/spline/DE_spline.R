#!/usr/bin/env Rscript
# =============================================================================
# DE_spline.R
#
# Differential expression across prenatal brain development, treating age
# (wpc) as a continuous variable. Each gene gets a smooth, non-linear
# trajectory (natural cubic spline), and we test whether that trajectory is
# flat. Age-varying genes are then clustered by the shape of their trajectory.
#
# Steps
#   1. Load raw counts + sample metadata; match columns to metadata rows.
#   2. Sanity checks (sample order, wpc in sample names, tissue vs age).
#   3. Design: intercept + ns(wpc, df = 3) [+ sex].
#   4. edgeR: filterByExpr -> TMM -> estimateDisp -> glmQLFit (robust).
#   5. Joint QL F-test on the spline coefficients = "changes with age".
#   6. Predict each gene's trajectory on a fine wpc grid; cluster the
#      significant genes (k-means on z-scored trajectories).
#   7. Save tables and figures.
#
# Inputs
#   quantification/Output/gene_counts_matrix_unique.csv   genes x samples, raw counts
#   analysis/data/sample_metadata.tsv                    sample_name, wpc, sex, ...
#
# Outputs
#   analysis/results/spline_DE/spline_DE_all_genes.csv
#   analysis/results/spline_DE/spline_DE_clusters.csv
#   analysis/results/spline_DE/qc_plots.pdf
#   Results-presentation/Figures/spline_trajectory_clusters.png
#   Results-presentation/Figures/spline_top_genes.png
#
# Run from the RStudio project (paths are built with here()).
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(splines)
  library(here)
  library(tidyverse)
})


# Settings ----
spline_df  <- 3      # flexibility of each gene's trajectory
adjust_sex <- FALSE   # include sex as a covariate
fdr_cutoff <- 0.05
n_clusters <- 6      # check the elbow plot in qc_plots.pdf and adjust
set.seed(1)

count_file <- here("quantification", "Output", "gene_counts_matrix_unique.csv")
meta_file  <- here("analysis", "data", "sample_metadata.tsv")
result_dir <- here("Results-presentation", "analysis-output")
figure_dir <- here("Results-presentation", "Figures")


# 1. Load data ----
dd   <- read_csv(count_file, show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

# Genes as row names; drop the "ERR..._" prefix from sample names
counts <- dd %>%
  column_to_rownames(var = names(dd)[1]) %>%
  as.matrix()
colnames(counts) <- sub("^[^_]*_", "", colnames(counts))


# 2. Match samples and sanity checks ----
missing <- setdiff(meta$sample_name, colnames(counts))
if (length(missing) > 0) {
  stop("In metadata but not in count matrix: ", paste(missing, collapse = ", "))
}
counts <- counts[, meta$sample_name]

# 3. Design matrix ----
age_basis <- ns(meta$wpc, df = spline_df)
colnames(age_basis) <- paste0("age_ns", seq_len(spline_df))

design <- cbind(Intercept = 1, age_basis)

rownames(design) <- meta$sample_name
spline_coefs <- 1 + seq_len(spline_df)   # columns tested for an age effect



# 4. edgeR: filter, normalise, fit ----
y <- DGEList(counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)                     # calcNormFactors() in edgeR < 4.0
message(sum(keep), " of ", length(keep), " genes kept after filtering")

y   <- estimateDisp(y, design, robust = TRUE)

fit <- glmQLFit(y, design, robust = TRUE)

pdf(file.path(result_dir, "qc_plots.pdf"), width = 7, height = 6)
plotMDS(y, pch = 16, col = hcl.colors(nlevels(factor(meta$wpc)), "Dark 3")[factor(meta$wpc)],
        main = "MDS (TMM-normalised logCPM)")
plotBCV(y, main = "Biological coefficient of variation")
plotQLDisp(fit, main = "QL dispersion")


# 5. Test: does expression change with age (any shape)? ----
res <- glmQLFTest(fit, coef = spline_coefs)
print(summary(decideTests(res, p.value = fdr_cutoff)))

tt <- topTags(res, n = Inf)$table %>%
  rownames_to_column("gene_id") %>%
  select(gene_id, logCPM, F, PValue, FDR)   # spline logFCs are not interpretable alone

write_csv(tt, file.path(result_dir, "spline_DE_all_genes.csv"))

sig_genes <- tt$gene_id[tt$FDR < fdr_cutoff]
message(length(sig_genes), " genes change with age at FDR < ", fdr_cutoff)


# 6. Trajectories on a fine wpc grid ----
# edgeR coefficients are natural-log expression relative to library size,
# so logCPM = (design %*% beta + log(1e6)) / log(2).
# Covariates other than age (e.g. sex) are set to their reference level;
# this shifts curves up/down only, not their shape.
wpc_grid   <- seq(min(meta$wpc), max(meta$wpc), length.out = 100)
grid_basis <- predict(age_basis, newx = wpc_grid)

grid_design <- matrix(0, nrow = length(wpc_grid), ncol = ncol(design),
                      dimnames = list(NULL, colnames(design)))
grid_design[, "Intercept"]  <- 1
grid_design[, spline_coefs] <- grid_basis

grid_logcpm <- (fit$coefficients %*% t(grid_design) + log(1e6)) / log(2)
colnames(grid_logcpm) <- wpc_grid

# Cluster significant genes by trajectory shape
if (length(sig_genes) <= n_clusters) {
  stop("Too few significant genes (", length(sig_genes), ") to cluster.")
}
z <- t(scale(t(grid_logcpm[sig_genes, , drop = FALSE])))

# Elbow plot to help choose n_clusters
wss <- map_dbl(2:12, \(k) kmeans(z, centers = k, nstart = 25, iter.max = 50)$tot.withinss)
plot(2:12, wss, type = "b", pch = 16, xlab = "Number of clusters (k)",
     ylab = "Total within-cluster SS", main = "Elbow plot for trajectory clustering")
dev.off()

km <- kmeans(z, centers = n_clusters, nstart = 50, iter.max = 50)

# Order clusters by when they peak (early -> late) for readable plots
peak_wpc <- apply(km$centers, 1, \(v) wpc_grid[which.max(v)])
cluster_order <- order(peak_wpc)
cluster_id <- match(km$cluster, cluster_order)
names(cluster_id) <- names(km$cluster)

cluster_tab <- tibble(gene_id = names(cluster_id), cluster = cluster_id) %>%
  left_join(tt, by = "gene_id") %>%
  arrange(cluster, FDR)
write_csv(cluster_tab, file.path(result_dir, "spline_DE_clusters.csv"))


# 7. Figures ----

# 7a. Trajectory clusters
traj_long <- as_tibble(z, rownames = "gene_id") %>%
  mutate(cluster = cluster_id[gene_id]) %>%
  pivot_longer(-c(gene_id, cluster), names_to = "wpc", values_to = "z") %>%
  mutate(wpc = as.numeric(wpc))

cluster_labels <- cluster_tab %>%
  dplyr::count(cluster) %>%
  mutate(label = paste0("Cluster ", cluster, " (n = ", n, ")"))
traj_long <- traj_long %>%
  mutate(cluster = factor(cluster, levels = cluster_labels$cluster,
                          labels = cluster_labels$label))

cluster_means <- traj_long %>%
  group_by(cluster, wpc) %>%
  summarise(z = mean(z), .groups = "drop")

p_clusters <- ggplot(traj_long, aes(wpc, z)) +
  geom_line(aes(group = gene_id), colour = "grey60", alpha = 0.05, linewidth = 0.3) +
  geom_line(data = cluster_means, colour = "#B2182B", linewidth = 1.2) +
  geom_rug(data = distinct(meta, wpc), aes(x = wpc), inherit.aes = FALSE,
           sides = "b", colour = "grey30") +
  facet_wrap(~ cluster) +
  labs(title = "Expression trajectories of age-varying genes",
       subtitle = paste0(length(sig_genes), " genes, FDR < ", fdr_cutoff,
                         "; spline df = ", spline_df),
       x = "Age (wpc)", y = "Scaled expression (z-score)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(figure_dir, "spline_trajectory_clusters.png"), p_clusters,
       width = 10, height = 6.5, dpi = 300, bg = "white")

# 7b. Top 9 genes: observed logCPM with fitted trajectory
top9 <- head(tt$gene_id, 9)
obs_logcpm <- cpm(y, log = TRUE, prior.count = 2)

obs_long <- as_tibble(obs_logcpm[top9, ], rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample_name", values_to = "logCPM") %>%
  left_join(select(meta, sample_name, wpc, sex), by = "sample_name")

fit_long <- as_tibble(grid_logcpm[top9, ], rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "wpc", values_to = "logCPM") %>%
  mutate(wpc = as.numeric(wpc))

gene_levels <- top9
s <- ggplot(obs_long, aes(wpc, logCPM)) +
  geom_point(aes(shape = sex), colour = "grey30", size = 2, alpha = 0.8) +
  geom_line(data = fit_long, colour = "#2166AC", linewidth = 1) +
  facet_wrap(~ factor(gene_id, levels = gene_levels), scales = "free_y") +
  labs(title = "Top age-varying genes",
       x = "Age (wpc)", y = "log2 CPM") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(figure_dir, "spline_top_genes.png"), p_top,
       width = 10, height = 8, dpi = 300, bg = "white")


message("Done. Results in ", result_dir, "; figures in ", figure_dir)
sessionInfo()
