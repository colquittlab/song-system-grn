## Topic-number metrics for the glut cells (hybrid labels): one small multiple per metric.
##   Rscript grn/ra-arco-hvc-nc_glut_hybrid/pycisTopic/plot_topic_metrics.R [out_dir]
## Reads models_grid_metrics.csv (pycistopic_glut_grid.py) and models_20topic_original_metrics.csv (the 20-topic model of
## models.pkl, trained at 40 threads; the grid's 20-topic model used 5, and mallet's parallel sampling is not deterministic
## across thread counts, so the second marker shows the run-to-run spread at K = 20).
suppressMessages({library(here); library(tidyverse)})
source(here::here("config/figure_theme.R"))

d <- here::here("grn/ra-arco-hvc-nc_glut_hybrid/pycisTopic")
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else
  "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_glut_hybrid/pycisTopic/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

metric_levels <- c("Arun_2010", "Cao_Juan_2009", "Mimno_2011", "loglikelihood")
panel_labels <- c(Arun_2010 = "Arun 2010\nlower is better",
                  Cao_Juan_2009 = "Cao-Juan 2009\nlower is better",
                  Mimno_2011 = "Mimno 2011\nhigher is better",
                  loglikelihood = "Log-likelihood (x1e8)\nhigher is better")
tidy <- function(f, run) read_csv(file.path(d, f), show_col_types = FALSE) %>%
  pivot_longer(-n_topic, names_to = "metric", values_to = "value") %>% mutate(run = run)
dat <- bind_rows(tidy("models_grid_metrics.csv", "Grid run (5 threads)"),
                 tidy("models_20topic_original_metrics.csv", "Original 20-topic model (40 threads)")) %>%
  mutate(value = if_else(metric == "loglikelihood", value / 1e8, value),
         metric = factor(metric, metric_levels, panel_labels[metric_levels]),
         run = factor(run, c("Grid run (5 threads)", "Original 20-topic model (40 threads)")))
grid <- filter(dat, run == "Grid run (5 threads)")
orig <- filter(dat, run != "Grid run (5 threads)")

p <- ggplot(grid, aes(n_topic, value)) +
  geom_line(colour = FIG_PAL[1], linewidth = 0.4) +
  geom_point(aes(colour = run, shape = run), size = 1.3) +
  geom_point(data = orig, aes(colour = run, shape = run), size = 1.6, stroke = 0.5) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  scale_colour_manual(NULL, values = c(FIG_PAL[1], FIG_PAL[2])) +
  scale_shape_manual(NULL, values = c(16, 1)) +
  scale_x_continuous("Number of topics", breaks = c(2, 10, 20, 30, 40), expand = fig_expand(0.04)) +
  scale_y_continuous(NULL, expand = fig_expand(0.06)) +
  theme_fig() +
  theme(strip.text = element_text(colour = FIG_INK_PRIMARY, size = FIG_PT_AXIS_TITLE, hjust = 0, face = "bold"),
        panel.spacing = unit(6, "pt"), legend.position = "bottom", legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = FIG_PT_AXIS_TEXT), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-2, 0, 0, 0))
fig_save(p, file.path(out_dir, "topic_metrics"), width = 3.4, height = 3.2)
cat("wrote", file.path(out_dir, "topic_metrics.{pdf,png}"), "\n")
