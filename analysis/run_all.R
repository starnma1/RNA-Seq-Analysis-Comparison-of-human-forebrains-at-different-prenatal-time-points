#!/usr/bin/env Rscript
# =============================================================================
# run_all.R
#
# Runs the whole R analysis, from sample metadata to GO enrichment, in order.
# Each script runs in its own fresh R process (like a restarted session), so
# no script can depend on objects or packages left behind by another one, or
# by whatever is already loaded in your RStudio session.
#
# Usage
#   In RStudio (from the project): source("analysis/run_all.R")
#   From a terminal, in the project root:  Rscript analysis/run_all.R
#
# Inputs
#   raw/samples.tsv, raw/ena/read_run.tsv, quantification/Output/*
#   (everything produced on Euler)
#
# Outputs
#   Everything the individual scripts write: analysis/data/,
#   Results-presentation/Figures/, Results-presentation/Tables/,
#   Results-presentation/analysis-output/
# =============================================================================

library(here)

# Settings ----
# Start from a later step to skip what is already up to date, e.g.
# run_from <- 8 reruns only the trajectory clustering and the enrichment.
run_from <- 1

steps <- c(
  "exploratory/01_samples.R",           #  1. sample metadata + overview
  "exploratory/02_multi_vs_unique.R",   #  2. featureCounts read assignment QC
  "Data-pre-processing.R",              #  3. protein-coding count matrix
  "exploratory/03_PCA.R",               #  4. PCA
  "exploratory/04_sex_check.R",         #  5. XIST vs RPS4Y1
  "exploratory/05_heatmap.R",           #  6. top 50 variable genes
  "spline/DE_spline.R",                 #  7. spline DE fit (slowest step)
  "spline/DE_spline_interpretation.R",  #  8. trajectory clusters + figures
  "GSEA/2_cluster_gsea.R"               #  9. GO enrichment by direction
)

# Run ----
# Start each script with the same R you are using now. "R" (not "Rscript")
# is used because it is always next to the running R, on every install.
r_bin <- file.path(R.home("bin"), "R")
if (!file.exists(r_bin)) stop("Can't find the R executable at ", r_bin, call. = FALSE)

# Run from the project root so here() inside each script finds the project
old_wd <- setwd(here())

for (i in seq(run_from, length(steps))) {
  script <- file.path("analysis", steps[i])
  message("\n[", i, "/", length(steps), "] ", steps[i])
  start <- Sys.time()

  # pdf(NULL): plots a script doesn't save itself are discarded instead of
  # ending up in an Rplots.pdf in the project root
  # --no-restore: don't load the project's .RData into the fresh session
  status <- system2(r_bin, c("--no-echo", "--no-restore", "--no-save", "-e", shQuote(
    sprintf("pdf(NULL); source('%s', echo = FALSE)", script)
  )))

  if (status != 0) {
    setwd(old_wd)
    stop("Step ", i, " (", steps[i], ") failed; see the error above. ",
         "Fix it, then set run_from <- ", i, " to continue from here.",
         call. = FALSE)
  }

  message("    done in ", format(round(Sys.time() - start, 1)))
}

setwd(old_wd)

message("\nAll steps finished.")
sessionInfo()
