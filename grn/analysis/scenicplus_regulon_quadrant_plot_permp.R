#!/usr/bin/env Rscript
# Same strip plot, point area scaled by -log10 of the permutation p-value of the UR share (scenicplus_regulon_quadrant_permutation.R: regulon-sized random sets
# from the pooled regulon genes / regions). Percentage of each regulon's genes (RNA log2FC) and regions (accessibility log2FC) in the upper-right quadrant of the RA vs C1H-1 / PVALB-1 vs LAMP5
# scatters (log2FC > 0 in both), one dot per regulon, MAFB highlighted with its rank (1 = highest percentage). Table: scenicplus_regulon_quadrant_fractions.R.
#
#   Rscript scenicplus_regulon_quadrant_plot_permp.R
suppressMessages({library(tidyverse); library(ggrepel); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
LGE <- identical(Sys.getenv("PV_COMPARATOR", "LAMP5"), "LGE1"); SUF <- if (LGE) "_lge1" else ""     # PV_COMPARATOR=LGE1: PVALB-1 vs GABA-LGE-1 instead of vs LAMP5
PLOT_DIR <- if (LGE) file.path(OUT, MAIN, "regulon_scatter_plots", "PVALB1_vs_LGE1") else file.path(OUT, MAIN, "regulon_scatter_plots")
R <- read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_fractions_", MAIN, SUF, ".csv")), show_col_types = FALSE) %>%
  mutate(pct = 100 * frac_UR, kind = factor(ifelse(kind == "gene_expression", "genes (RNA log2FC)", "regions (accessibility log2FC)"),
                                            levels = c("genes (RNA log2FC)", "regions (accessibility log2FC)"))) %>%
  left_join(read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_permutation_", MAIN, SUF, ".csv")), show_col_types = FALSE) %>%
              transmute(TF, kind = ifelse(kind == "gene_expression", "genes (RNA log2FC)", "regions (accessibility log2FC)"), p_frac, nlp = -log10(p_frac)), by = c("TF", "kind")) %>%
  mutate(kind = factor(kind, levels = c("genes (RNA log2FC)", "regions (accessibility log2FC)"))) %>%
  group_by(kind) %>% mutate(rank = rank(-pct, ties.method = "min"), n_reg = n()) %>% ungroup()
set.seed(2026)
R <- R %>% mutate(xi = as.numeric(kind), xj = ifelse(TF == "MAFB", xi, xi + runif(n(), -0.18, 0.18)))   # one fixed jitter for every dot
M <- R %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB\n%.1f%%, rank %d of %d\npermutation p = %.2g", pct, rank, n_reg, p_frac))
print(as.data.frame(M %>% select(kind, pct, rank, n_reg)))
# regulons with a higher percentage than MAFB in the same version are labeled
A <- R %>% left_join(M %>% select(kind, pct_mafb = pct), by = "kind") %>% filter(TF != "MAFB", pct > pct_mafb)
print(as.data.frame(A %>% arrange(kind, desc(pct)) %>% select(kind, TF, pct, n) %>% mutate(pct = round(pct, 1))))
p <- ggplot(R, aes(xj, pct)) +
  geom_point(data = R %>% filter(TF != "MAFB"), aes(size = nlp), color = FIG_INK_MUTED, alpha = 0.75, stroke = 0) +
  geom_point(data = A, aes(size = nlp), color = FIG_INK_PRIMARY, stroke = 0) +
  geom_text_repel(data = A, aes(label = TF), nudge_x = -0.4, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, min.segment.length = 0,
                  hjust = 1, direction = "y", seed = 1, box.padding = 0.3, force = 1.5, max.time = 5) +
  geom_point(data = M, aes(size = nlp), color = FIG_PAL[2], stroke = 0) +
  geom_text_repel(data = M, aes(label = lab), nudge_x = 0.4, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, min.segment.length = 0,
                  lineheight = 0.95, hjust = 0, direction = "y") +
  scale_size_continuous(name = expression(-log[10]~p[perm]), range = c(0.3, 2.4), breaks = c(1, 2, 3)) +
  scale_x_continuous(breaks = 1:2, labels = sub(" \\(", "\n(", levels(R$kind)), limits = c(0.1, 3.3), expand = c(0, 0)) +
  labs(x = NULL, y = "% in the upper-right quadrant\n(log2FC > 0 in both contrasts)") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TITLE), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
stem <- file.path(PLOT_DIR, "quadrant_percentages_by_regulon_perm_p")
fig_save(p, stem, width = 4.2, height = 3.6)
ggsave(paste0(stem, ".png"), p, width = 4.2, height = 3.6, dpi = 600, bg = "transparent")
ggsave(paste0(stem, ".svg"), p, width = 4.2, height = 3.6, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
