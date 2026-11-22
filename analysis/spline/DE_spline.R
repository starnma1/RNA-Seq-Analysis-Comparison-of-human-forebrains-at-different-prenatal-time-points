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
#   1. Load the protein-coding count matrix + sample metadata; match
#      columns to metadata rows.
#   2. Check every sample is present and put columns in metadata order.
#   3. Design: intercept + ns(wpc, df = 3)
#   4. edgeR: filterByExpr -> TMM -> estimateDisp -> glmQLFit (robust).
#   5. Joint QL F-test on the spline coefficients = "changes with age".
#   6. Save the fit (for DE_spline_interpretation.R) and the results table.
#
# Inputs
#   analysis/data/processed_count_matrix.csv   protein-coding genes x samples, raw
#                                              counts (from Data-pre-processing.R)
#   analysis/data/sample_metadata.tsv          run, sample_name, wpc, sex, ...
#
# Outputs (Results-presentation/analysis-output/)
#   spline_fit.rds             fitted model + everything the later scripts need
#   spline_DE_all_genes.csv    spline test result for every tested gene
#   qc_plots.pdf               MDS, BCV and QL dispersion plots
#
# Clustering and figures are done in DE_spline_interpretation.R.
#
# Run from the RStudio project (paths are built with here()).
# =============================================================================

  library(edgeR)
  library(splines)
  library(here)
  library(tidyverse)


# Settings ----
spline_df  <- 3      # flexibility of each gene's trajectory
adjust_sex <- FALSE   # include sex as a covariate
fdr_cutoff <- 0.05
n_clusters <- 6      # check the elbow plot in qc_plots.pdf and adjust

count_file <- here("analysis", "data", "processed_count_matrix.csv")   # protein-coding only
meta_file  <- here("analysis", "data", "sample_metadata.tsv")
result_dir <- here("Results-presentation", "analysis-output")
figure_dir <- here("Results-presentation", "Figures")

###############################################################################

# 1. Load data ----
dd   <- read_csv(count_file, show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

# Genes as row names; columns are run IDs (ERR...)
counts <- dd %>%
  column_to_rownames(var = names(dd)[1]) %>%
  as.matrix()


# 2. Match samples ----
# Put columns in metadata order, then rename them to sample names
# (DE_spline_interpretation.R joins counts and metadata by sample_name).
stopifnot(setequal(colnames(counts), meta$run))
counts <- counts[, meta$run]
colnames(counts) <- meta$sample_name

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

y   <- estimateDisp(y, design, robust = TRUE)

fit <- glmQLFit(y, design, robust = TRUE)

# Quality checks
pdf(file.path(result_dir, "qc_plots.pdf"), width = 7, height = 6)

pal <- hcl.colors(100, "Viridis", rev = TRUE)
wpc_idx <- round((meta$wpc - min(meta$wpc)) / diff(range(meta$wpc)) * 99) + 1

mds <- plotMDS(
  y,
  pch = 16,
  cex = 1.4,
  col = pal[wpc_idx],
  main = "MDS (TMM-normalised logCPM)"
)

leg_weeks <- c(4, 8, 12, 16, 20)
leg_idx   <- round((leg_weeks - min(meta$wpc)) / diff(range(meta$wpc)) * 99) + 1
legend("topright", legend = paste(leg_weeks, "wpc"),
       pch = 16, pt.cex = 1.3, col = pal[leg_idx], bty = "n", title = "Age")

plotBCV(y, main = "Biological coefficient of variation")
plotQLDisp(fit, main = "QL dispersion")

dev.off()   # <-- Closes and finalizes qc_plots.pdf


# 5. Test: does expression change with age (any shape)? ----
res <- glmQLFTest(fit, coef = spline_coefs)
print(summary(decideTests(res, p.value = fdr_cutoff)))

tt <- topTags(res, n = Inf)$table %>%
  rownames_to_column("gene_id") %>%
  select(gene_id, logCPM, F, PValue, FDR)   # spline logFCs are not interpretable alone


sig_genes <- tt$gene_id[tt$FDR < fdr_cutoff]

write_csv(tt, file.path(result_dir, "spline_DE_all_genes.csv"))


# 6. Save fit for the interpretation script ----
saveRDS(
  list(fit          = fit,
       y            = y,             # needed for observed logCPM in the top-genes plot
       age_basis    = age_basis,
       design       = design,
       meta         = meta,
       tt           = tt,
       sig_genes    = sig_genes,
       spline_coefs = spline_coefs,
       fdr_cutoff   = fdr_cutoff,
       spline_df    = spline_df),
  file.path(result_dir, "spline_fit.rds")
)

sessionInfo()
