#!/usr/bin/env Rscript
# Flipped-axis version of the compact UR-share strip plot (scenicplus_regulon_quadrant_plot_permp_compact.R): percentage on x, genes and regions as two
# horizontal strips packed close together. Point area scales with -log10 permutation p of the UR share (scenicplus_regulon_quadrant_permutation.R),
# MAFB orange, regulons above MAFB labeled in a column to the right of the strips. 2.6 x 1.6 in.
#
#   Rscript scenicplus_regulon_quadrant_plot_permp_flipped.R
suppressMessages({library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
KIND <- c(gene_expression = "genes", region_accessibility = "regions")
YPOS <- c(genes = 1.55, regions = 1)      # strip centers (genes on top); the gap between strips is deliberately small
R <- read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_fractions_", MAIN, ".csv")), show_col_types = FALSE) %>%
  select(TF, kind, n, frac_UR) %>%
  left_join(read_csv(here::here(paste0("grn/analysis/scenicplus_regulon_quadrant_permutation_", MAIN, ".csv")), show_col_types = FALSE) %>% select(TF, kind, p_frac), by = c("TF", "kind")) %>%
  mutate(pct = 100 * frac_UR, nlp = -log10(p_frac), kind = factor(KIND[kind], levels = KIND), yi = YPOS[as.character(kind)]) %>%
  group_by(kind) %>% mutate(rank = rank(-pct, ties.method = "min")) %>% ungroup()
set.seed(2026)
R <- R %>% mutate(yj = ifelse(TF == "MAFB", yi, yi + runif(n(), -0.16, 0.16)))
M <- R %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB #%d, %.1f%%", rank, pct))
A <- R %>% left_join(M %>% select(kind, pct_mafb = pct), by = "kind") %>% filter(TF != "MAFB", pct > pct_mafb)
# labeled dots (regulons above MAFB, then MAFB) are staggered in y in share order, highest on top, so no leader line crosses another labeled dot
ST <- bind_rows(A, M %>% select(any_of(names(A)))) %>% group_by(kind) %>% arrange(desc(pct), .by_group = TRUE) %>% mutate(yj = yi + 0.17 - 0.055 * (row_number() - 1)) %>% ungroup()
A <- ST %>% filter(TF != "MAFB"); M <- M %>% select(-yj) %>% left_join(ST %>% filter(TF == "MAFB") %>% select(kind, yj), by = "kind")
# label column right of the strips: stacked in share order, minimum vertical spacing
LX <- 56; GAP <- 0.095
LAB <- bind_rows(A %>% mutate(lab = TF), M) %>% arrange(kind, desc(pct)) %>% group_by(kind) %>%
  mutate(ly = { y <- numeric(n()); for (k in seq_len(n())) y[k] <- if (k == 1) yj[k] + (if (kind[1] == "genes") 0.06 else -0.05) else min(yj[k], y[k - 1] - GAP); y }) %>% ungroup()
p <- ggplot(R, aes(pct, yj)) +
  geom_point(data = R %>% filter(TF != "MAFB", !paste(TF, kind) %in% paste(A$TF, A$kind)), aes(size = nlp), color = FIG_INK_MUTED, alpha = 0.75, stroke = 0) +
  geom_segment(data = LAB, aes(x = LX - 0.8, xend = pct, y = ly, yend = yj), linewidth = 0.2, color = FIG_INK_PRIMARY) +
  geom_point(data = A, aes(size = nlp), color = FIG_INK_PRIMARY, stroke = 0) +
  geom_point(data = M, aes(size = nlp), color = FIG_PAL[2], stroke = 0) +
  geom_text(data = LAB, aes(x = LX, y = ly, label = lab), hjust = 0, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY) +
  scale_size_continuous(name = expression(-log[10]~p), range = c(0.2, 1.6), breaks = c(1, 2, 3)) +
  scale_x_continuous(breaks = seq(10, 50, 10), limits = c(5, 77), expand = c(0, 0)) +
  scale_y_continuous(breaks = unname(YPOS), labels = names(YPOS), limits = c(0.82, 1.86), expand = c(0, 0)) +
  labs(x = "% in upper-right quadrant", y = NULL) + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.12, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0), legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
stem <- file.path(OUT, MAIN, "regulon_scatter_plots", "quadrant_percentages_by_regulon_perm_p_flipped")
W <- 2.6; H <- 1.6
fig_save(p, stem, width = W, height = H)
ggsave(paste0(stem, ".png"), p, width = W, height = H, dpi = 600, bg = "transparent")
ggsave(paste0(stem, ".svg"), p, width = W, height = H, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
