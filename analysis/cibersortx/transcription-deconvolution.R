#!/usr/bin/env Rscript
# =============================================================================
# plot_deconvolution.R
#
# Plots BrainDeconvShiny results (neuron / astrocyte fractions) against age.
#
# Inputs (analysis/data/)
#   braindeconv_CIBERSORT.csv   downloaded from the app ("CIBERSORT" button)
#   braindeconv_dtangle.csv     downloaded from the app ("DTANGLE" button), optional
#   sample_metadata.tsv         run, sample_name, sex, wpc, ...
#
# Outputs (Results-presentation/Figures/)
#   deconvolution_vs_wpc.png        fraction per cell type vs wpc, with loess trend
#   deconvolution_per_sample.png    stacked bars, samples ordered by age
# =============================================================================

library(here)
library(tidyverse)

data_dir   <- here("analysis", "data")
figure_dir <- here("Results-presentation", "Figures")

result_files <- c(
  CIBERSORT = file.path(data_dir, "CIBERSORT.csv"),
  dtangle   = file.path(data_dir, "dtangle.csv")
)
result_files <- result_files[file.exists(result_files)]
if (!length(result_files)) stop("No result files found in ", data_dir)

cell_types <- c("Neurons", "Astrocytes", "Oligodendrocytes",
                "Microglia", "Endothelia", "OPCs", "Excitatory", "Inhibitory")

meta <- read_tsv(file.path(data_dir, "sample_metadata.tsv"), show_col_types = FALSE)

# ---- Read results: first column = sample, keep cell-type columns only ------
read_result <- function(file, method) {
  read_csv(file, show_col_types = FALSE) %>%
    rename(run = 1) %>%
    select(run, any_of(cell_types)) %>%
    pivot_longer(-run, names_to = "cell_type", values_to = "fraction") %>%
    mutate(method = method)
}

fractions <- imap_dfr(result_files, read_result) %>%
  inner_join(meta, by = "run")

stopifnot(nrow(fractions) > 0)   # fails if sample IDs don't match metadata$run
message(n_distinct(fractions$run), " samples, methods: ",
        paste(unique(fractions$method), collapse = ", "))

# ---- 1. Fraction vs age -----------------------------------------------------
p1 <- ggplot(fractions, aes(wpc, fraction, colour = cell_type)) +
  geom_point(aes(shape = sex), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ method) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Age (wpc)", y = "Estimated fraction", colour = "Cell type",
       shape = "Sex",
       title = "Estimated cell-type composition across development",
       subtitle = "BrainDeconvShiny, adult cortical signature") +
  theme_bw(base_size = 12)

ggsave(file.path(figure_dir, "deconvolution_vs_wpc.png"), p1,
       width = 4 + 3.5 * length(result_files), height = 4.5, dpi = 300)

# ---- 2. Stacked bars per sample, ordered by age -----------------------------
sample_order <- meta %>% arrange(wpc, run) %>% pull(run)

p2 <- fractions %>%
  mutate(run = factor(run, levels = sample_order)) %>%
  ggplot(aes(run, fraction, fill = cell_type)) +
  geom_col(width = 0.9) +
  facet_wrap(~ method, ncol = 1) +
  scale_x_discrete(labels = setNames(paste0(meta$wpc, " wpc"), meta$run)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Estimated fraction", fill = "Cell type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

ggsave(file.path(figure_dir, "deconvolution_per_sample.png"), p2,
       width = 10, height = 3 + 2.5 * length(result_files), dpi = 300)

# ---- Quick check: does the neuron fraction rise with age? -------------------
fractions %>%
  filter(cell_type == "Neurons") %>%
  group_by(method) %>%
  summarise(spearman_rho = cor(wpc, fraction, method = "spearman"),
            p_value      = cor.test(wpc, fraction, method = "spearman",
                                    exact = FALSE)$p.value)