#!/usr/bin/env Rscript
# =============================================================================
# DE_spline_interpretation.R
#
# Interpretation of the spline DE fit from DE_spline.R. Rebuilds each gene's
# fitted trajectory on a fine wpc grid, clusters the age-varying genes by the
# shape of their trajectory, and makes the figures.
#
# No model fitting happens here, so clustering and plots can be changed and
# rerun quickly without refitting edgeR.
#
# Steps
#   1. Load the saved fit (spline_fit.rds) from DE_spline.R.
#   2. Predict each gene's trajectory (log2 CPM) on a fine wpc grid.
#   3. Choose k: elbow plot + preview of cluster shapes for several k.
#   4. Final k-means clustering (z-scored trajectories), clusters ordered
#      by peak age (early -> late).
#   5. Figures: trajectory clusters, top 9 genes with fitted curves.
#
# Inputs
#   Results-presentation/analysis-output/spline_fit.rds   saved by DE_spline.R
#
# Outputs
#   Results-presentation/analysis-output/spline_DE_clusters.csv
#   Results-presentation/Figures/spline_elbow.png
#   Results-presentation/Figures/spline_cluster_preview.png
#   Results-presentation/Figures/spline_trajectory_clusters.png
#   Results-presentation/Figures/spline_top_genes.png
#
# Run from the RStudio project (paths are built with here()).
# =============================================================================

library(edgeR)
library(splines)
library(here)
library(tidyverse)


# Settings ----
n_clusters <- 6      # set after checking the elbow and cluster preview plots
k_try      <- 4:8    # k values shown in the cluster preview
n_grid     <- 100    # number of ages on the prediction grid
seed       <- 1

result_dir <- here("Results-presentation", "analysis-output")
figure_dir <- here("Results-presentation", "Figures")
fit_file   <- file.path(result_dir, "spline_fit.rds")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################

# 1. Load the saved fit ----
sp <- readRDS(fit_file)

fit          <- sp$fit
y            <- sp$y
age_basis    <- sp$age_basis
design       <- sp$design
meta         <- sp$meta
tt           <- sp$tt
sig_genes    <- sp$sig_genes
spline_coefs <- sp$spline_coefs
fdr_cutoff   <- sp$fdr_cutoff
spline_df    <- sp$spline_df

message(length(sig_genes), " age-varying genes (FDR < ", fdr_cutoff, ")")

if (length(sig_genes) <= max(c(n_clusters, k_try))) {
  stop("Too few significant genes (", length(sig_genes), ") to cluster.")
}


# 2. Trajectories on a fine wpc grid ----
# edgeR coefficients are natural-log expression relative to library size,
# so logCPM = (design %*% beta + log(1e6)) / log(2).
# Covariates other than age (e.g. sex) are set to their reference level;
# this shifts curves up/down only, not their shape.
wpc_grid   <- seq(min(meta$wpc), max(meta$wpc), length.out = n_grid)
grid_basis <- predict(age_basis, newx = wpc_grid)   # same knots as the fit

grid_design <- matrix(0, nrow = length(wpc_grid), ncol = ncol(design),
                      dimnames = list(NULL, colnames(design)))
grid_design[, "Intercept"]  <- 1
grid_design[, spline_coefs] <- grid_basis

grid_logcpm <- (fit$coefficients %*% t(grid_design) + log(1e6)) / log(2)
colnames(grid_logcpm) <- wpc_grid

# Scale each gene's curve (mean 0, SD 1) so clustering uses shape only
z <- t(scale(t(grid_logcpm[sig_genes, , drop = FALSE])))


# 3. Choosing the number of clusters ----
# 3a. Elbow plot
set.seed(seed)
k_elbow <- 2:12
wss <- sapply(k_elbow, function(k) {
  kmeans(z, centers = k, nstart = 25, iter.max = 50)$tot.withinss
})

p_elbow <- ggplot(tibble(k = k_elbow, wss = wss), aes(k, wss)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = k_elbow) +
  labs(title = "Elbow plot for trajectory clustering",
       x = "Number of clusters (k)", y = "Total within-cluster SS") +
  theme_bw()
print(p_elbow)

ggsave(file.path(figure_dir, "spline_elbow.png"), p_elbow,
       width = 6, height = 4.5, dpi = 300, bg = "white")

# 3b. Cluster shapes for several k
set.seed(seed)
cluster_preview <- data.frame()   # empty, filled in the loop

for (k in k_try) {
  km_k <- kmeans(z, centers = k, nstart = 25, iter.max = 50)
  
  result_k <- data.frame(
    gene_id = rownames(z),
    cluster = km_k$cluster,
    k = k
  )
  
  cluster_preview <- rbind(cluster_preview, result_k)   # add this k's rows
}

z_long <- as_tibble(z, rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "wpc", values_to = "z") %>%
  mutate(wpc = as.numeric(wpc))

preview_means <- cluster_preview %>%
  left_join(z_long, by = "gene_id", relationship = "many-to-many") %>%
  group_by(k, cluster, wpc) %>%
  summarise(mean_z  = mean(z),
            lo      = quantile(z, 0.25),
            hi      = quantile(z, 0.75),
            n_genes = n(),
            .groups = "drop")

# Number clusters by peak age within each k, so rows read early -> late
preview_order <- preview_means %>%
  group_by(k, cluster) %>%
  summarise(peak = wpc[which.max(mean_z)], .groups = "drop") %>%
  arrange(k, peak) %>%
  group_by(k) %>%
  mutate(panel = row_number()) %>%
  ungroup()

preview_means <- preview_means %>%
  left_join(preview_order, by = c("k", "cluster"))

preview_sizes <- preview_means %>%
  distinct(k, panel, n_genes)

p_preview <- ggplot(preview_means, aes(wpc, mean_z)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
  geom_line(colour = "steelblue4", linewidth = 0.8) +
  geom_text(data = preview_sizes, aes(label = paste0("n = ", n_genes)),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3) +
  facet_grid(paste("k =", k) ~ panel) +
  labs(title = "Cluster shapes for different numbers of clusters",
       subtitle = "Line = cluster mean; band = middle 50% of genes",
       x = "Age (wpc)", y = "Scaled expression (z)") +
  theme_bw()
print(p_preview)

ggsave(file.path(figure_dir, "spline_cluster_preview.png"), p_preview,
       width = 14, height = 10, dpi = 300, bg = "white")


# 4. Final clustering ----
set.seed(seed)
km <- kmeans(z, centers = n_clusters, nstart = 50, iter.max = 50)

# Order clusters by when they peak (early -> late) for readable plots
peak_wpc      <- apply(km$centers, 1, function(v) wpc_grid[which.max(v)])
cluster_order <- order(peak_wpc)
cluster_id    <- match(km$cluster, cluster_order)
names(cluster_id) <- names(km$cluster)

cluster_tab <- tibble(gene_id = names(cluster_id), cluster = cluster_id) %>%
  left_join(tt, by = "gene_id") %>%
  arrange(cluster, FDR)

write_csv(cluster_tab, file.path(result_dir, "spline_DE_clusters.csv"))
print(table(cluster = cluster_tab$cluster))


# 5. Figures ----

# 5a. Trajectory clusters
cluster_labels <- cluster_tab %>%
  dplyr::count(cluster) %>%
  mutate(label = paste0("Cluster ", cluster, " (n = ", n, ")"))

traj_long <- z_long %>%
  mutate(cluster = cluster_id[gene_id],
         cluster = factor(cluster, levels = cluster_labels$cluster,
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
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_clusters)

ggsave(file.path(figure_dir, "spline_trajectory_clusters.png"), p_clusters,
       width = 10, height = 6.5, dpi = 300, bg = "white")

# 5b. Top 9 genes: observed logCPM with fitted trajectory
top9 <- head(tt$gene_id, 9)
obs_logcpm <- cpm(y, log = TRUE, prior.count = 2)

obs_long <- as_tibble(obs_logcpm[top9, ], rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample_name", values_to = "logCPM") %>%
  left_join(select(meta, sample_name, wpc, sex), by = "sample_name")

fit_long <- as_tibble(grid_logcpm[top9, ], rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "wpc", values_to = "logCPM") %>%
  mutate(wpc = as.numeric(wpc))

p_top <- ggplot(obs_long, aes(wpc, logCPM)) +
  geom_point(aes(shape = sex), colour = "grey30", size = 2, alpha = 0.8) +
  geom_line(data = fit_long, colour = "#2166AC", linewidth = 1) +
  facet_wrap(~ factor(gene_id, levels = top9), scales = "free_y") +
  labs(title = "Top age-varying genes",
       x = "Age (wpc)", y = "log2 CPM") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_top)

ggsave(file.path(figure_dir, "spline_top_genes.png"), p_top,
       width = 10, height = 8, dpi = 300, bg = "white")


message("Done. Results in ", result_dir, "; figures in ", figure_dir)
sessionInfo()