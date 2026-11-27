#!/usr/bin/env Rscript
# marker_check.R — is the "endothelia" fraction really endothelium, or progenitors?

library(here)
library(tidyverse)
library(edgeR)

data_dir   <- here("analysis", "data")
figure_dir <- here("Results-presentation", "Figures")

meta   <- read_tsv(file.path(data_dir, "sample_metadata.tsv"), show_col_types = FALSE)
genes  <- read_tsv(file.path(data_dir, "gene_annotation.tsv"), show_col_types = FALSE)
counts <- as.matrix(read.csv(file.path(data_dir, "protein_coding_counts.csv"),
                             row.names = 1, check.names = FALSE))[, meta$run]

markers <- tribble(
  ~gene,    ~group,
  "PECAM1", "Endothelial", "CDH5", "Endothelial", "CLDN5", "Endothelial",
  "SOX2",   "Progenitor",  "PAX6", "Progenitor",  "NES",   "Progenitor",
  "MKI67",  "Progenitor",
  "STMN2",  "Neuronal",    "SNAP25","Neuronal"
)

# TMM-normalised logCPM
y <- normLibSizes(DGEList(counts))
logcpm <- cpm(y, log = TRUE, prior.count = 2)
rownames(logcpm) <- genes$gene_name[match(rownames(logcpm), genes$gene_id)]

expr <- logcpm[intersect(markers$gene, rownames(logcpm)), , drop = FALSE] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "run", values_to = "logCPM") %>%
  left_join(markers, by = "gene") %>%
  left_join(meta, by = "run")

p <- ggplot(expr, aes(wpc, logCPM, colour = gene)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_smooth(method = "loess", span = 1, se = FALSE) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "Age (wpc)", y = "log2 CPM (TMM)", colour = "Gene") +
  theme_bw(base_size = 12)
ggsave(file.path(figure_dir, "marker_check.png"), p, width = 11, height = 4, dpi = 300)

# Which markers track the "Endothelia" fraction?
endo <- read_csv(file.path(data_dir, "CIBERSORT.csv"), show_col_types = FALSE) %>%
  rename(run = 1) %>% select(run, Endothelia)

expr %>%
  inner_join(endo, by = "run") %>%
  group_by(group, gene) %>%
  summarise(rho_with_endothelia = cor(logCPM, Endothelia, method = "spearman"),
            .groups = "drop") %>%
  arrange(desc(rho_with_endothelia)) %>%
  print()
