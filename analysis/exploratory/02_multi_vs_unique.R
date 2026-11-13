#!/usr/bin/env Rscript
# =============================================================================
# compare_unique_multi.R
#
# Summarises what happened to the reads in featureCounts (unique-read counting):
# how many were assigned to genes, how many were set aside as multimappers,
# and how many fell outside genes or were ambiguous. This is the justification
# for using the unique-read count matrix in all downstream analyses.
#
# Input
#   quantification/Output/featurecounts_summary_unique.tsv
#
# Outputs
#   analysis/data/read_fate_unique.tsv
#   Results-presentation/Tables/read_fate_unique.png
#       Median and range (across all samples) of the % of reads in each
#       featureCounts category.
# =============================================================================

library(here)
library(tidyverse)
library(gridExtra)
library(grid)

data_dir  <- here("quantification", "Output")
out_dir  <- here("analysis", "data")
table_dir <- here("Results-presentation", "Tables")

# ---- Read the featureCounts summary (rows = read categories, columns = samples) ----
fc_unique <- read_tsv(file.path(data_dir, "featurecounts_summary_unique.tsv"))

# ---- % of reads per category, summarised across samples ----------------------
read_fate <- fc_unique %>%
  pivot_longer(-Status, names_to = "sample", values_to = "reads") %>%
  group_by(sample) %>%
  mutate(pct = 100 * reads / sum(reads)) %>%
  group_by(Status) %>%
  summarise(median_pct = round(median(pct), 1),
            min_pct    = round(min(pct), 1),
            max_pct    = round(max(pct), 1),
            .groups = "drop") %>%
  filter(median_pct > 0) %>%                  # drop categories that are always 0
  arrange(desc(median_pct))

# ---- Save as an image, sized to fit the table --------------------------------------
tbl <- tableGrob(
  read_fate %>% rename(`Read category` = Status,
                       `Median %` = median_pct,
                       `Min %` = min_pct,
                       `Max %` = max_pct),
  rows = NULL, theme = ttheme_minimal(base_size = 11)
)

ggsave(file.path(table_dir, "read_fate_unique.png"), tbl,
       width  = convertWidth(sum(tbl$widths), "in", valueOnly = TRUE) + 0.4,
       height = convertHeight(sum(tbl$heights), "in", valueOnly = TRUE) + 0.4,
       dpi = 300, bg = "white")