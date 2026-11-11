#!/usr/bin/env Rscript
# =============================================================================
# sample_metadata.R
#
# Describes the 32 prenatal brain samples (E-MTAB-6814) used in this project.
#
# Inputs
#   raw/samples.tsv           run ID, sample name, BioSample, wpc, sex
#   raw/ena/read_run.tsv      ENA run report (library details, read counts)
#
# Outputs
#   1. analysis/data/sample_metadata.tsv
#        Metadata table for all downstream analyses (one row per sample):
#        run, sample_name, sex, wpc, library_layout, library_selection,
#        read_count, read_length.
#        Use it as colData for DESeq2 (match rows to count-matrix columns by run).
#   2. analysis/data/sample_summary.tsv  +  Results-presentation/Tables/sample_summary.png
#        Summary table: number of samples (total, female, male) and median
#        sequencing depth per developmental age (wpc).
#   3. Results-presentation/Figures/sample_overview.png
#        Figure: every sample by age (x) and sequencing depth (y), coloured by sex.
#
# Notes
#   All samples are prenatal 'brain' in E-MTAB-6814 = forebrain + midbrain,
#   whether labelled 'Forebrain' or 'Brain' in the sample name.
#   Run from the RStudio project (paths are built with here()).
# =============================================================================

library(here)
library(tidyverse)
library(gridExtra) 
library(grid)

# Output folders (created if they don't exist yet)
data_dir   <- here("analysis", "data")
table_dir  <- here("Results-presentation", "Tables")
figure_dir <- here("Results-presentation", "Figures")

# ---- Read inputs -----------------------------------------------------------------
samples  <- read_tsv(here("raw", "samples.tsv"), show_col_types = FALSE)
read_run <- read_tsv(here("raw", "ena", "read_run.tsv"), show_col_types = FALSE)

# ---- 1. Metadata table -----------------------------------------------------------
metadata <- samples %>%
  left_join(read_run, by = c("run" = "run_accession")) %>%
  transmute(
    run,
    sample_name,
    sex,
    wpc,
    library_layout,
    library_selection,
    read_count,
    read_length = base_count / read_count    # average read length (bp)
  ) %>%
  arrange(wpc, sex, run)

write_tsv(metadata, file.path(data_dir, "sample_metadata.tsv"))

# ---- 2. Summary table ------------------------------------------------------------
summary_table <- metadata %>%
  group_by(wpc) %>%
  summarise(
    n_samples      = n(),
    n_female       = sum(sex == "female"),
    n_male         = sum(sex == "male"),
    median_reads_M = round(median(read_count) / 1e6, 1),
    .groups = "drop"
  )

write_tsv(summary_table, file.path(data_dir, "sample_summary.tsv"))

# Save the table as an image, sized to fit the table exactly
tbl <- tableGrob(
  summary_table %>%
    rename(wpc = wpc, Samples = n_samples, Female = n_female, Male = n_male,
           `Median reads (M)` = median_reads_M),
  rows = NULL, theme = ttheme_minimal(base_size = 11)
)
ggsave(file.path(table_dir, "sample_summary.png"), tbl,
       width  = convertWidth(sum(tbl$widths), "in", valueOnly = TRUE) + 0.4,
       height = convertHeight(sum(tbl$heights), "in", valueOnly = TRUE) + 0.4,
       dpi = 300)

# ---- 3. Figure -------------------------------------------------------------------
p <- ggplot(metadata) +
  
  geom_point(aes(x = wpc, y = read_count / 1e6, colour = sex),
             size = 3,
             position = position_jitter(width = 0.15, height = 0)) +
  
  scale_colour_manual(values = c(female = "#eb6834", male = "#2a78d6")) +
  scale_x_continuous(breaks = sort(unique(metadata$wpc))) +
  scale_y_continuous(limits = c(0, NA)) +
  
  labs(title = "Samples by developmental age",
       x = "Weeks post-conception (wpc)",
       y = "Sequenced reads (millions)",
       colour = "Sex") +
  
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Save the file
ggsave(file.path(figure_dir, "sample_overview.png"), p,
       width = 8, height = 4.5, dpi = 300)