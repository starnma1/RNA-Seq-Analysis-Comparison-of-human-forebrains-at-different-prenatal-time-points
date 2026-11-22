# =============================================================================
# Over-representation analysis (GO BP) of genes that change with age
#
# Takes the edgeR spline fit and splits the significant genes by the direction
# of their fitted expression change with wpc, then runs GO enrichment on each
# group. Two ways of defining direction:
#   1. Endpoints:   predicted expression at the oldest vs the youngest wpc
#   2. Trajectory:  Spearman correlation of the full fitted curve with wpc,
#                   which also separates out transient (non-monotonic) genes
#
# Input:  spline_fit.rds (from the spline fitting script)
# Output: dot plots in Figures/, enrichment tables in analysis-output/
# =============================================================================

# Load libraries ----
library(here)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)   # change if not human

# Settings ----
n_grid        <- 100    # number of ages on the prediction grid
trend_cutoff  <- 0.8    # |Spearman rho| needed to call a gene increasing/decreasing
go_ontology   <- "BP"   # "BP", "MF" or "CC"
show_category <- 10     # GO terms shown per group in the dot plots
fig_width     <- 8      # inches
fig_height    <- 9      # inches

group_order <- c("Decreasing", "Increasing", "Transient")   # column order in plots

fit_file   <- here("Results-presentation", "analysis-output", "spline_fit.rds")
result_dir <- here("Results-presentation", "analysis-output")
figure_dir <- here("Results-presentation", "Figures")

# Load data needed ----
sp <- readRDS(fit_file)

fit          <- sp$fit
design       <- sp$design
meta         <- sp$meta
age_basis    <- sp$age_basis
spline_coefs <- sp$spline_coefs
sig_genes    <- sp$sig_genes

# Helper functions ----

# Drop the Ensembl version suffix (ENSG00000141510.17 -> ENSG00000141510)
strip <- function(x) sub("\\..*$", "", x)

# Fitted log2 expression for the given genes at the given ages
# (rows = genes, columns = ages). Builds a design matrix with one row per age,
# using the same spline basis (same knots) as the original fit.
predict_expression <- function(ages, genes) {
  X <- matrix(0, nrow = length(ages), ncol = ncol(design),
              dimnames = list(NULL, colnames(design)))
  X[, "Intercept"]  <- 1
  X[, spline_coefs] <- predict(age_basis, newx = ages)
  
  out <- (fit$coefficients[genes, , drop = FALSE] %*% t(X)) / log(2)
  colnames(out) <- ages
  out
}

# GO over-representation for each gene set against all tested genes,
# with groups always in group_order. Saves the table and the dot plot.
run_ora <- function(gene_sets, name) {
  gene_sets <- gene_sets[intersect(group_order, names(gene_sets))]
  
  ora <- compareCluster(gene_sets, fun = "enrichGO",
                        OrgDb    = org.Hs.eg.db,
                        keyType  = "ENSEMBL",
                        ont      = go_ontology,
                        universe = strip(rownames(fit)))
  
  ora@compareClusterResult$Cluster <- factor(ora@compareClusterResult$Cluster,
                                             levels = names(gene_sets))
  
  write.csv(as.data.frame(ora),
            file.path(result_dir, paste0(name, ".csv")), row.names = FALSE)
  
  p <- dotplot(ora, showCategory = show_category) +
    theme(axis.text.y = element_text(size = 7))
  
  ggsave(file.path(figure_dir, paste0(name, ".png")), p,
         width = fig_width, height = fig_height, dpi = 300)
  
  p
}

# Method 1: direction from the endpoints (youngest vs oldest wpc) ----
ends   <- range(meta$wpc)
pred   <- predict_expression(ends, sig_genes)
change <- pred[, 2] - pred[, 1]   # log2 change, oldest minus youngest

gene_sets_ends <- list(
  Decreasing = strip(sig_genes[change < 0]),
  Increasing = strip(sig_genes[change > 0])
)
lengths(gene_sets_ends)

first_vs_last_ora <- run_ora(gene_sets_ends, "first_vs_last_ora")
first_vs_last_ora

# Method 2: direction from the full fitted trajectory ----
wpc_grid <- seq(min(meta$wpc), max(meta$wpc), length.out = n_grid)
traj     <- predict_expression(wpc_grid, sig_genes)

# Spearman rho of each curve with age: +1 = only rises, -1 = only falls
trend <- apply(traj, 1, function(g) cor(g, wpc_grid, method = "spearman"))

direction <- case_when(
  trend >=  trend_cutoff ~ "Increasing",
  trend <= -trend_cutoff ~ "Decreasing",
  TRUE                   ~ "Transient"
)
names(direction) <- rownames(traj)
table(direction)

gene_sets_traj <- split(strip(names(direction)), direction)

trajectory_ora <- run_ora(gene_sets_traj, "trajectory_inc_vs_dec_ora")
trajectory_ora