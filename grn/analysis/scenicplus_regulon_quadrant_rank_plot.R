#!/usr/bin/env Rscript
# Percentage of each regulon's genes (RNA log2FC) and regions (accessibility log2FC) in the upper-right quadrant, plotted against its rank (1 = highest
# percentage), genes and regions in separate panels. MAFB highlighted; regulons above MAFB labeled; point area scales with -log10 permutation p of the UR
# share (scenicplus_regulon_quadrant_permutation.R). Table: scenicplus_regulon_quadrant_fractions.R.
#
#   Rscript scenicplus_regulon_quadrant_rank_plot.R
suppressMessages({library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
LGE <- identical(Sys.getenv("PV_COMPARATOR", "LAMP5"), "LGE1"); SUF <- if (LGE) "_lge1" else ""     # PV_COMPARATOR=LGE1: PVALB-1 vs GABA-LGE-1 instead of vs LAMP5
PLOT_DIR <- if (LGE) file.path(OUT, MAIN, "regulon_scatter_plots", "PVALB1_vs_LGE1") else file.path(OUT, MAIN, "regulon_scatter_plots")
KIND <- c(gene_expression = "genes (RNA log2FC)", region_accessibility = "regions (accessibility log2FC)")
R <- read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_fractions_", MAIN, SUF, ".csv")), show_col_types = FALSE) %>%
  select(TF, kind, n, frac_UR) %>%
  left_join(read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_permutation_", MAIN, SUF, ".csv")), show_col_types = FALSE) %>% select(TF, kind, p_frac), by = c("TF", "kind")) %>%
  mutate(pct = 100 * frac_UR, nlp = -log10(p_frac), kind = factor(KIND[kind], levels = KIND)) %>%
  group_by(kind) %>% arrange(desc(pct), .by_group = TRUE) %>% mutate(rank = rank(-pct, ties.method = "first")) %>% ungroup()
M <- R %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB: #%d, %.1f%%", rank, pct))
A <- R %>% filter(TF != "MAFB") %>% left_join(M %>% select(kind, rank_mafb = rank), by = "kind") %>% filter(rank < rank_mafb)
# labels stacked in rank order to the right of the top-ranked points (no crossing leader lines)
LAB <- bind_rows(A %>% mutate(lab = TF), M) %>% arrange(kind, rank) %>% group_by(kind) %>%
  mutate(ly = { y <- numeric(n()); for (k in seq_len(n())) y[k] <- if (k == 1) pct[k] else min(pct[k], y[k - 1] - 3.6); y }, lx = 26) %>% ungroup()
p <- ggplot(R, aes(rank, pct)) +
  geom_point(data = R %>% filter(TF != "MAFB", !paste(TF, kind) %in% paste(A$TF, A$kind)), aes(size = nlp), color = FIG_INK_MUTED, alpha = 0.75, stroke = 0) +
  geom_point(data = A, aes(size = nlp), color = FIG_INK_PRIMARY, stroke = 0) +
  geom_point(data = M, aes(size = nlp), color = FIG_PAL[2], stroke = 0) +
  geom_segment(data = LAB, aes(x = lx - 0.8, xend = rank, y = ly, yend = pct), linewidth = 0.2, color = FIG_INK_PRIMARY) +
  geom_text(data = LAB, aes(x = lx, y = ly, label = lab), hjust = 0, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY) +
  facet_wrap(~kind, nrow = 1) +
  scale_size_continuous(name = expression(-log[10]~p), range = c(0.2, 1.6), breaks = c(1, 2, 3)) +
  scale_x_continuous(expand = fig_expand(), breaks = c(1, 20, 40, 60, 79)) + scale_y_continuous(expand = fig_expand()) +
  labs(x = "Rank of regulon (1 = highest %)", y = "% in upper-right quadrant") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(), legend.key.size = unit(0.15, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT),
        strip.text = element_text(size = FIG_PT_AXIS_TITLE, family = FIG_FONT), strip.background = element_blank(), panel.spacing = unit(0.08, "in"))
stem <- file.path(PLOT_DIR, "quadrant_percentage_by_rank")
W <- 3.6; H <- 1.6
fig_save(p, stem, width = W, height = H)
ggsave(paste0(stem, ".png"), p, width = W, height = H, dpi = 600, bg = "transparent")
ggsave(paste0(stem, ".svg"), p, width = W, height = H, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
