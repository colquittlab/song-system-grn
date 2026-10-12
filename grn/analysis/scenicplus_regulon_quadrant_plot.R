#!/usr/bin/env Rscript
# Percentage of each regulon's genes (RNA log2FC) and regions (accessibility log2FC) in the upper-right quadrant of the RA vs C1H-1 / PVALB-1 vs LAMP5
# scatters (log2FC > 0 in both), one dot per regulon, MAFB highlighted with its rank (1 = highest percentage). Table: scenicplus_regulon_quadrant_fractions.R.
#
#   Rscript scenicplus_regulon_quadrant_plot.R
suppressMessages({library(tidyverse); library(ggrepel); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
R <- read_csv(here::here("grn/analysis/scenicplus_regulon_quadrant_fractions_config37.csv"), show_col_types = FALSE) %>%
  mutate(pct = 100 * frac_UR, kind = factor(ifelse(kind == "gene_expression", "genes (RNA log2FC)", "regions (accessibility log2FC)"),
                                            levels = c("genes (RNA log2FC)", "regions (accessibility log2FC)"))) %>%
  group_by(kind) %>% mutate(rank = rank(-pct, ties.method = "min"), n_reg = n()) %>% ungroup()
set.seed(2026)
R <- R %>% mutate(xi = as.numeric(kind), xj = ifelse(TF == "MAFB", xi, xi + runif(n(), -0.18, 0.18)))   # one fixed jitter for every dot
M <- R %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB\n%.1f%%, rank %d of %d", pct, rank, n_reg))
print(as.data.frame(M %>% select(kind, pct, rank, n_reg)))
# regulons with a higher percentage than MAFB in the same version are labeled
A <- R %>% left_join(M %>% select(kind, pct_mafb = pct), by = "kind") %>% filter(TF != "MAFB", pct > pct_mafb)
print(as.data.frame(A %>% arrange(kind, desc(pct)) %>% select(kind, TF, pct, n) %>% mutate(pct = round(pct, 1))))
p <- ggplot(R, aes(xj, pct)) +
  geom_point(data = R %>% filter(TF != "MAFB"), size = 0.9, color = FIG_INK_MUTED, alpha = 0.75, stroke = 0) +
  geom_point(data = A, size = 0.9, color = FIG_INK_PRIMARY, stroke = 0) +
  geom_text_repel(data = A, aes(label = TF), nudge_x = -0.4, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, min.segment.length = 0,
                  hjust = 1, direction = "y", seed = 1, box.padding = 0.3, force = 1.5, max.time = 5) +
  geom_point(data = M, size = 2.2, color = FIG_PAL[2], stroke = 0) +
  geom_text_repel(data = M, aes(label = lab), nudge_x = 0.4, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, min.segment.length = 0,
                  lineheight = 0.95, hjust = 0, direction = "y") +
  scale_x_continuous(breaks = 1:2, labels = levels(R$kind), limits = c(0.1, 3.3), expand = c(0, 0)) +
  labs(x = NULL, y = "% in the upper-right quadrant\n(log2FC > 0 in both contrasts)") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank())
stem <- file.path(OUT, MAIN, "regulon_scatter_plots", "quadrant_percentages_by_regulon")
fig_save(p, stem, width = 3.3, height = 3.6)
ggsave(paste0(stem, ".png"), p, width = 3.3, height = 3.6, dpi = 600, bg = "transparent")
ggsave(paste0(stem, ".svg"), p, width = 3.3, height = 3.6, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
