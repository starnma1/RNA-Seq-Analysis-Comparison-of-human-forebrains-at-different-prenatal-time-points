#!/usr/bin/env Rscript
# =============================================================================
# PCA.R
#
# Exploratory PCA of the prenatal brain samples, to check whether
# developmental age is the main source of variation in the data.
#
# Steps
#   1. Variance-stabilizing transformation (DESeq2::vst) of the filtered,
#      protein-coding count matrix (for exploration only; DE uses raw counts).
#   2. Select the 500 most variable genes.
#   3. PCA with samples as observations.
#   4. Plot PC1 vs PC2: samples coloured by age (wpc), plus the mean position
#      of 5 age groups (4–5, 7–8, 9–11, 12–13, 16–19 wpc) connected in
#      developmental order.
#
# Inputs (analysis/data/)
#   processed_count_matrix.csv   genes x samples, raw counts (from Data-pre-processing.R)
#   gene_annotation.tsv          gene names and biotypes
#   sample_metadata.tsv          run, wpc, sex, ... (from sample_metadata.R)
#
# Output
#   Results-presentation/Figures/pca_age_groups.png
#
# Run from the RStudio project (paths are built with here()).
# =============================================================================

library(DESeq2)
library(here)
library(tidyverse)

count_dir  <- here("analysis", "data")
figure_dir <- here("Results-presentation", "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


# Load in the Data ----
counts  <- read.csv(file.path(count_dir, "processed_count_matrix.csv"),
                    row.names = 1, check.names = FALSE) %>% as.matrix()

genes   <- read.csv(file.path(count_dir, "gene_annotation.tsv"), sep = "\t")
samples <- read.csv(file.path(count_dir, "sample_metadata.tsv"), sep = "\t")


# Variance-stabilizing transformation
vst_mat <- vst(counts)

# Top 500 most variable genes
gene_variance <- apply(vst_mat, 1, var)
top500 <- names(sort(gene_variance, decreasing = TRUE))[1:500]

# PCA (samples as rows)
pca <- prcomp(t(vst_mat[top500, ]))
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1) 

# PCA scores + metadata, with 5 age groups
df <- data.frame(pca$x[, 1:2], wpc = samples$wpc, sex = samples$sex) %>%
  mutate(age_group = cut(wpc, breaks = c(0, 5, 8, 11, 13, Inf),
                         labels = c("4–5 wpc", "7–8 wpc", "9–11 wpc", "12–13 wpc", "16–19 wpc")))

# Mean PC1/PC2 per age group
group_means <- df %>%
  group_by(age_group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = n())

# Plot
p <- ggplot(df, aes(PC1, PC2)) +
  
  geom_point(aes(colour = wpc, shape = sex), size = 3.25, alpha = 0.8) +
  geom_point(data = group_means, shape = 21, fill = "white", colour = "black", size = 5) +
  
  geom_text(data = group_means, aes(label = age_group),
            vjust = -1.4, size = 3.5, fontface = "bold", alpha = .3) +
  
  labs(title = "PCA of Prenatal Forebrain samples",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p
ggsave(file.path(figure_dir, "pca_age_groups.png"), p,
       width = 8, height = 6, dpi = 300, bg = "white")
