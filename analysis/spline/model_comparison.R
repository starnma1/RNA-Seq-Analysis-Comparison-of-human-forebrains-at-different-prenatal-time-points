#!/usr/bin/env Rscript
# =============================================================================
# model_comparison.R
#
# Is a straight line in log(wpc) enough, or do genes need a curved trajectory?
#
#   1. Orthogonal cubic polynomial in log2(wpc): coefficient 1 = straight-line
#      trend, coefficients 2-3 = curvature. A QL F-test on the curvature terms
#      finds genes that a line in log(wpc) cannot describe.
#   2. Sensitivity analysis: repeat the curvature test without the youngest
#      samples (wpc < min_wpc_sensitivity), to see how much of the curvature
#      depends on the earliest developmental window.
#   3. Plot the top curvature genes: observed logCPM, straight line in
#      log(wpc), and the natural spline used in DE_spline.R.
#
# Inputs (analysis/data/)
#   protein_coding_counts.csv, sample_metadata.tsv, gene_annotation.tsv
#
# Outputs
#   Results-presentation/analysis-output/model_comparison_summary.tsv
#   Results-presentation/analysis-output/model_comparison_curvature_genes.csv
#                                       (includes the sensitivity-analysis FDR)
#   Results-presentation/Figures/model_comparison_top_curvature_genes.png
# =============================================================================

library(edgeR)
library(splines)
library(here)
library(tidyverse)

data_dir   <- here("analysis", "data")
result_dir <- here("Results-presentation", "analysis-output")
figure_dir <- here("Results-presentation", "Figures")
fdr_cutoff <- 0.05
n_show     <- 6      # top curvature genes to plot
spline_df  <- 3      # same as DE_spline.R
min_wpc_sensitivity <- 7   # sensitivity analysis drops samples younger than this

# ---- Load ----------------------------------------------------------------------
meta   <- read_tsv(file.path(data_dir, "sample_metadata.tsv"), show_col_types = FALSE)
genes  <- read_tsv(file.path(data_dir, "gene_annotation.tsv"), show_col_types = FALSE)
counts <- as.matrix(read.csv(file.path(data_dir, "protein_coding_counts.csv"),
                             row.names = 1, check.names = FALSE))
counts <- counts[, meta$run]

log_wpc <- log2(meta$wpc)

# ---- Designs -------------------------------------------------------------------
# Polynomial in log2(wpc): orthogonal, so "linear" is the straight-line part
# and "quad"/"cubic" are curvature on top of it
design_poly <- model.matrix(~ poly(log_wpc, 3))
colnames(design_poly) <- c("Intercept", "linear", "quad", "cubic")

# Straight line in log2(wpc) and natural spline in wpc, for plotting
design_line <- cbind(Intercept = 1, log_wpc = log_wpc)
age_basis   <- ns(meta$wpc, df = spline_df)
design_spl  <- cbind(Intercept = 1, age_basis)

# ---- Filter and normalise once, so every model uses the same genes ----------------
y <- DGEList(counts)
y <- y[filterByExpr(y, design = design_poly), , keep.lib.sizes = FALSE]
y <- normLibSizes(y)

fit_model <- function(design, dge = y) {
  dge <- estimateDisp(dge, design, robust = TRUE)
  glmQLFit(dge, design, robust = TRUE)
}

fit_poly <- fit_model(design_poly)
fit_line <- fit_model(design_line)
fit_spl  <- fit_model(design_spl)

# ---- Tests ---------------------------------------------------------------------
lin   <- glmQLFTest(fit_poly, coef = "linear")
curve <- glmQLFTest(fit_poly, coef = c("quad", "cubic"))

lin_fdr   <- p.adjust(lin$table$PValue,   "BH")
curve_fdr <- p.adjust(curve$table$PValue, "BH")

# ---- Sensitivity analysis: without the youngest samples -------------------------
# The 4-5 wpc samples sit at the edge of the age range and change fastest, so
# they have high leverage on the curvature test. Refit without them (own
# filtering and normalisation) and repeat the test.
keep_s      <- meta$wpc >= min_wpc_sensitivity
log_wpc_s   <- log2(meta$wpc[keep_s])
design_poly_s <- model.matrix(~ poly(log_wpc_s, 3))
colnames(design_poly_s) <- c("Intercept", "linear", "quad", "cubic")

y_s <- DGEList(counts[, keep_s])
y_s <- y_s[filterByExpr(y_s, design = design_poly_s), , keep.lib.sizes = FALSE]
y_s <- normLibSizes(y_s)

fit_poly_s  <- fit_model(design_poly_s, dge = y_s)
curve_s     <- glmQLFTest(fit_poly_s, coef = c("quad", "cubic"))
curve_fdr_s <- setNames(p.adjust(curve_s$table$PValue, "BH"), rownames(y_s))

curved_all <- rownames(y)[curve_fdr < fdr_cutoff]
curved_s   <- names(curve_fdr_s)[curve_fdr_s < fdr_cutoff]

message(sprintf("Sensitivity (wpc >= %d, %d of %d samples): %d genes curved (vs %d with all samples), %d in both",
                min_wpc_sensitivity, sum(keep_s), length(keep_s),
                length(curved_s), length(curved_all),
                length(intersect(curved_all, curved_s))))

summary_tab <- tibble(
  category = c("Genes tested",
               "Linear trend in log(wpc), up",
               "Linear trend in log(wpc), down",
               "Significant curvature (quad + cubic)",
               "Curvature but no linear trend",
               "Linear trend only (no curvature)",
               sprintf("Sensitivity (wpc >= %d): samples used", min_wpc_sensitivity),
               sprintf("Sensitivity (wpc >= %d): genes tested", min_wpc_sensitivity),
               sprintf("Sensitivity (wpc >= %d): significant curvature", min_wpc_sensitivity),
               "Curved with all samples AND without the youngest"),
  n_genes = c(nrow(y),
              sum(lin_fdr < fdr_cutoff & lin$table$logFC > 0),
              sum(lin_fdr < fdr_cutoff & lin$table$logFC < 0),
              sum(curve_fdr < fdr_cutoff),
              sum(curve_fdr < fdr_cutoff & lin_fdr >= fdr_cutoff),
              sum(lin_fdr < fdr_cutoff & curve_fdr >= fdr_cutoff),
              sum(keep_s),
              nrow(y_s),
              length(curved_s),
              length(intersect(curved_all, curved_s)))
)
print(summary_tab)
write_tsv(summary_tab, file.path(result_dir, "model_comparison_summary.tsv"))

curvature_genes <- topTags(curve, n = Inf)$table %>%
  rownames_to_column("gene_id") %>%
  left_join(select(genes, gene_id, gene_name), by = "gene_id") %>%
  mutate(FDR_without_youngest = unname(curve_fdr_s[gene_id]),          # NA = not tested
         curved_without_youngest = !is.na(FDR_without_youngest) &
                                   FDR_without_youngest < fdr_cutoff) %>%
  select(gene_id, gene_name, logCPM, F, PValue, FDR,
         FDR_without_youngest, curved_without_youngest)
write_csv(curvature_genes, file.path(result_dir, "model_comparison_curvature_genes.csv"))

# ---- Fitted curves on a fine age grid ----------------------------------------------
wpc_grid <- seq(min(meta$wpc), max(meta$wpc), length.out = 200)

# edgeR coefficients are on the natural-log scale -> convert to log2 CPM
to_logcpm <- function(fit, X, ids) {
  (fit$coefficients[ids, , drop = FALSE] %*% t(X) + log(1e6)) / log(2)
}

top_ids <- head(curvature_genes$gene_id, n_show)

X_line <- cbind(1, log2(wpc_grid))
X_spl  <- cbind(1, predict(age_basis, newx = wpc_grid))

curves <- bind_rows(
  as_tibble(to_logcpm(fit_line, X_line, top_ids), rownames = "gene_id") %>%
    mutate(model = "Straight line in log(wpc)"),
  as_tibble(to_logcpm(fit_spl, X_spl, top_ids), rownames = "gene_id") %>%
    mutate(model = "Natural spline (df = 3)")
) %>%
  pivot_longer(-c(gene_id, model), names_to = "i", values_to = "logCPM") %>%
  mutate(wpc = wpc_grid[as.integer(sub("^V", "", i))]) %>%
  select(-i)

# Observed TMM-normalised logCPM
observed <- cpm(y, log = TRUE, prior.count = 2)[top_ids, , drop = FALSE] %>%
  as_tibble(rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "run", values_to = "logCPM") %>%
  left_join(select(meta, run, wpc, sex), by = "run")

# Panel labels: gene name + curvature FDR, in order of significance
labels <- curvature_genes %>%
  filter(gene_id %in% top_ids) %>%
  mutate(label = sprintf("%s\ncurvature FDR = %.1e",
                         coalesce(gene_name, gene_id), FDR))
lab <- setNames(labels$label, labels$gene_id)

curves   <- mutate(curves,   gene = factor(lab[gene_id], levels = lab[top_ids]))
observed <- mutate(observed, gene = factor(lab[gene_id], levels = lab[top_ids]))

# ---- Plot ----------------------------------------------------------------------
p <- ggplot() +
  geom_point(data = observed, aes(wpc, logCPM, shape = sex),
             colour = "grey35", size = 2, alpha = 0.8) +
  geom_line(data = curves, aes(wpc, logCPM, colour = model, linetype = model),
            linewidth = 1) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("Straight line in log(wpc)" = "#d62728",
                                 "Natural spline (df = 3)"  = "#1f77b4")) +
  scale_linetype_manual(values = c("Straight line in log(wpc)" = "dashed",
                                   "Natural spline (df = 3)"  = "solid")) +
  labs(x = "Age (wpc)", y = "log2 CPM (TMM)", colour = NULL, linetype = NULL,
       shape = "Sex",
       title = "Genes a straight line in log(wpc) cannot describe",
       subtitle = sprintf(paste0("%d of %d genes show significant curvature (FDR < %.2f); ",
                                 "%d remain curved without samples < %d wpc"),
                          length(curved_all), nrow(y), fdr_cutoff,
                          length(curved_s), min_wpc_sensitivity)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(figure_dir, "model_comparison_top_curvature_genes.png"), p,
       width = 10, height = 7, dpi = 300)