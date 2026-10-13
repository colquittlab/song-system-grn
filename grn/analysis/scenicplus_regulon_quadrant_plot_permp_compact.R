#!/usr/bin/env Rscript
# Compact (2.6 x 1.6 in total) version of the same strip plot, point area scaled by -log10 of the permutation p-value of the UR share (scenicplus_regulon_quadrant_permutation.R: regulon-sized random sets
# from the pooled regulon genes / regions). Percentage of each regulon's genes (RNA log2FC) and regions (accessibility log2FC) in the upper-right quadrant of the RA vs C1H-1 / PVALB-1 vs LAMP5
# scatters (log2FC > 0 in both), one dot per regulon, MAFB highlighted with its rank (1 = highest percentage). Table: scenicplus_regulon_quadrant_fractions.R.
#
#   Rscript scenicplus_regulon_quadrant_plot_permp_compact.R
suppressMessages({library(tidyverse); library(ggrepel); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
R <- read_csv(here::here("grn/analysis/scenicplus_regulon_quadrant_fractions_config37.csv"), show_col_types = FALSE) %>%
  mutate(pct = 100 * frac_UR, kind = factor(ifelse(kind == "gene_expression", "genes (RNA log2FC)", "regions (accessibility log2FC)"),
                                            levels = c("genes (RNA log2FC)", "regions (accessibility log2FC)"))) %>%
  left_join(read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_permutation_", MAIN, ".csv")), show_col_types = FALSE) %>%
              transmute(TF, kind = ifelse(kind == "gene_expression", "genes (RNA log2FC)", "regions (accessibility log2FC)"), p_frac, nlp = -log10(p_frac)), by = c("TF", "kind")) %>%
  mutate(kind = factor(kind, levels = c("genes (RNA log2FC)", "regions (accessibility log2FC)"))) %>%
  group_by(kind) %>% mutate(rank = rank(-pct, ties.method = "min"), n_reg = n()) %>% ungroup()
set.seed(2026)
R <- R %>% mutate(xi = as.numeric(kind), xj = ifelse(TF == "MAFB", xi, xi + runif(n(), -0.18, 0.18)))   # one fixed jitter for every dot
M <- R %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB\n%.1f%%, rank %d", pct, rank))
print(as.data.frame(M %>% select(kind, pct, rank, n_reg)))
# regulons with a higher percentage than MAFB in the same version are labeled
A <- R %>% left_join(M %>% select(kind, pct_mafb = pct), by = "kind") %>% filter(TF != "MAFB", pct > pct_mafb)
# labeled dots are staggered in x (highest share leftmost) so no leader line crosses another labeled dot; their grey duplicates are dropped below
A <- A %>% group_by(kind) %>% arrange(desc(pct), .by_group = TRUE) %>% mutate(xj = xi - 0.16 + 0.05 * (row_number() - 1)) %>% ungroup()
print(as.data.frame(A %>% arrange(kind, desc(pct)) %>% select(kind, TF, pct, n) %>% mutate(pct = round(pct, 1))))
# label positions: gene labels (regulons above MAFB, and MAFB) share one column left of the genes strip, stacked top-down with a minimum spacing;
# the region labels sit beside (MAFB) or above (NR3C2) their dots. Leader lines run from the label to the dot.
GAP <- 3.4
LAB <- bind_rows(A %>% mutate(lab = TF), M %>% mutate(lab = sprintf("MAFB\n#%d, %.1f%%", rank, pct))) %>% arrange(kind, desc(pct))
g1 <- LAB %>% filter(kind == levels(R$kind)[1]); ly <- numeric(nrow(g1)); prev_gap <- 0
for (k in seq_len(nrow(g1))) { need <- if (k == 1) Inf else ly[k - 1] - ifelse(g1$TF[k] == "MAFB", GAP * 1.8, GAP); ly[k] <- min(g1$pct[k], need) }
g1 <- g1 %>% mutate(ly = ly, lx = 0.32, hj = 0, sx = lx + 0.065 * pmax(nchar(TF), ifelse(TF == "MAFB", 7, 0)) + 0.05, sy = ly)
g2 <- LAB %>% filter(kind == levels(R$kind)[2]) %>% mutate(ly = ifelse(TF == "MAFB", pct, pct + 5), lx = ifelse(TF == "MAFB", 2.3, 1.75), hj = ifelse(TF == "MAFB", 0, 0.5),
                                                          sx = ifelse(TF == "MAFB", 2.27, 1.75), sy = ifelse(TF == "MAFB", pct, ly - 1.6))
LAB <- bind_rows(g1, g2)
p <- ggplot(R, aes(xj, pct)) +
  geom_point(data = R %>% filter(TF != "MAFB", !paste(TF, kind) %in% paste(A$TF, A$kind)), aes(size = nlp), color = FIG_INK_MUTED, alpha = 0.75, stroke = 0) +
  geom_point(data = A, aes(size = nlp), color = FIG_INK_PRIMARY, stroke = 0) +
  geom_segment(data = LAB, aes(x = sx, xend = xj, y = sy, yend = pct), linewidth = 0.2, color = FIG_INK_PRIMARY) +
  geom_text(data = LAB, aes(x = lx, y = ly, label = lab, hjust = hj), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, lineheight = 0.9) +
  geom_point(data = M, aes(size = nlp), color = FIG_PAL[2], stroke = 0) +
  scale_size_continuous(name = expression(-log[10]~p), range = c(0.2, 1.6), breaks = c(1, 2, 3)) +
  scale_x_continuous(breaks = 1:2, labels = c("genes", "regions"), limits = c(0.3, 3.1), expand = c(0, 0)) +
  labs(x = NULL, y = "% in upper-right quadrant") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.15, "cm"), legend.margin = margin(0, 0, 0, 0), legend.position = "inside", legend.position.inside = c(0.93, 0.2), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
stem <- file.path(OUT, MAIN, "regulon_scatter_plots", "quadrant_percentages_by_regulon_perm_p_compact")
fig_save(p, stem, width = 2.6, height = 1.6)
ggsave(paste0(stem, ".png"), p, width = 2.6, height = 1.6, dpi = 600, bg = "transparent")
ggsave(paste0(stem, ".svg"), p, width = 2.6, height = 1.6, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
