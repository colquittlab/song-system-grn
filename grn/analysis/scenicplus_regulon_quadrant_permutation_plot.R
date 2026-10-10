#!/usr/bin/env Rscript
# Permutation test of the UR ratio (scenicplus_regulon_quadrant_permutation.R): observed UR ratio of each regulon against regulon size, with the permutation
# null (random sets of the same size from the pooled regulon genes / regions) as a band of mean +/- 1.96 SD. Point area scales with -log10 (empirical p of the ratio),
# MAFB highlighted, regulons with BH q < 0.05 labeled.
#
#   Rscript scenicplus_regulon_quadrant_permutation_plot.R
suppressMessages({library(tidyverse); library(ggrepel); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
R <- read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_permutation_", MAIN, ".csv")), show_col_types = FALSE) %>%
  mutate(kind = factor(ifelse(kind == "gene_expression", "genes (RNA log2FC)", "regions (accessibility log2FC)"), levels = c("genes (RNA log2FC)", "regions (accessibility log2FC)")))
L <- bind_rows(
  R %>% transmute(TF, kind, n, metric = "UR ratio", obs = ratio, nm = null_ratio_mean, nsd = null_ratio_sd, p = p_ratio, q = q_ratio),
  R %>% transmute(TF, kind, n, metric = "UR share", obs = frac_UR, nm = null_frac_mean, nsd = null_frac_sd, p = p_frac, q = q_frac)) %>%
  mutate(nlp = -log10(p), metric = factor(metric, c("UR ratio", "UR share")))
M <- L %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB\np = %.2g", p))
# labels: ratio row = BH q < 0.05; share row = regulons with a higher share than MAFB (the same set labeled in quadrant_percentages_by_regulon)
S <- L %>% left_join(M %>% select(kind, metric, obs_mafb = obs), by = c("kind", "metric")) %>% filter(TF != "MAFB", ifelse(metric == "UR ratio", q < 0.05, obs > obs_mafb))
p <- ggplot(L, aes(n, obs)) +
  geom_ribbon(data = L %>% arrange(n), aes(ymin = nm - 1.96 * nsd, ymax = nm + 1.96 * nsd), fill = FIG_GRID, alpha = 0.8) +
  geom_line(data = L %>% arrange(n), aes(y = nm), linewidth = 0.3, color = FIG_INK_MUTED) +
  geom_point(data = L %>% filter(TF != "MAFB"), aes(size = nlp), color = FIG_INK_MUTED, alpha = 0.7, stroke = 0) +
  geom_point(data = S, aes(size = nlp), color = FIG_INK_PRIMARY, stroke = 0) +
  geom_text_repel(data = S, aes(label = TF), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, min.segment.length = 0, box.padding = 0.4, seed = 1, max.overlaps = Inf) +
  geom_point(data = M, aes(size = nlp), color = FIG_PAL[2], stroke = 0) +
  geom_text_repel(data = M, aes(label = lab), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, min.segment.length = 0, lineheight = 0.95, box.padding = 0.5, seed = 1) +
  facet_grid(metric ~ kind, scales = "free", switch = "y") +
  scale_size_continuous(name = expression(-log[10]~p), range = c(0.4, 3.2), breaks = c(0, 1, 2, 3)) +
  scale_x_log10(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
  labs(x = "Regulon size (genes or regions)", y = "Observed (grey band: permutation null, mean ± 1.96 SD)") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TITLE), legend.text = element_text(size = FIG_PT_AXIS_TEXT),
        strip.text = element_text(size = FIG_PT_AXIS_TITLE, family = FIG_FONT), strip.background = element_blank(), strip.placement = "outside", panel.spacing = unit(0.2, "in"))
stem <- file.path(OUT, MAIN, "regulon_scatter_plots", "quadrant_permutation")
fig_save(p, stem, width = 6.2, height = 4.6)
ggsave(paste0(stem, ".svg"), p, width = 6.2, height = 4.6, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
