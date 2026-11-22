#!/usr/bin/env Rscript
# =============================================================================
# 04_sex_check.R
#
# Checks that each sample's annotated sex matches its expression of a
# female-specific gene (XIST) and a male-specific gene (RPS4Y1).
#
# Input
#   quantification/Output/gene_counts_annotated_unique.csv
#       gene_id, gene_name, gene_biotype, gene_length, then one column of raw
#       counts per sample (sample names end in Female/Male)
#
# Output
#   Scatter plot of XIST vs RPS4Y1 counts, coloured by annotated sex
#   (shown in the plot pane, not saved)
#
# Run from the RStudio project (paths are built with here()).
# =============================================================================

library(here)

# Load data ----
gene_counts_annotated_unique <- read.csv(
  here("quantification", "Output", "gene_counts_annotated_unique.csv"),
  check.names = FALSE
)

all_genes <- gene_counts_annotated_unique
xist_row <- all_genes[all_genes$gene_name == "XIST", ]
rps4y1_row <- all_genes[all_genes$gene_name == "RPS4Y1", ]

stopifnot(nrow(xist_row) == 1, nrow(rps4y1_row) == 1)
stopifnot(identical(names(xist_row), names(rps4y1_row)))

x <- unlist(xist_row[1, 5:ncol(xist_row)], use.names = FALSE)
y <- unlist(rps4y1_row[1, 5:ncol(rps4y1_row)], use.names = FALSE)

samples <- names(xist_row)[5:ncol(xist_row)]
colors <- ifelse(grepl("Female", samples, ignore.case = TRUE), "red",
                 ifelse(grepl("Male", samples, ignore.case = TRUE), "blue", "grey"))

plot(x + 1, y + 1, log = "xy", col = colors, pch = 19,
     xlab = "XIST counts + 1", ylab = "RPS4Y1 counts + 1")

legend("topright",
       legend = c(
         paste0("Female (n=", sum(colors == "red"), ")"),
         paste0("Male (n=", sum(colors == "blue"), ")")
       ),
       col = c("red", "blue"), pch = 19)

