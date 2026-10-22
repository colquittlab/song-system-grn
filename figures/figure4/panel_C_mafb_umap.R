#!/usr/bin/env Rscript
# Figure 4, panel C -- MAFB expression on the multiome UMAP.
#
# Each cell on the hybrid-labeled object's UMAP (reduction dims25nn30mindist0.3: Harmony-integrated PCA, 25 dimensions, 30 neighbours, min.dist 0.3;
# multiome/seurat/reduction_viz/combined_all_umap_hybrid.R). Cells without MAFB are light gray; cells with MAFB are drawn on top in increasing order of
# expression, colored by SCT log-normalized expression on a single-hue scale capped at the 99th percentile of the expressing cells (so a few extreme cells do
# not wash out the rest). Points are rasterized at 600 dpi. 2.2 x 2.55 in.
# (Earlier versions of the figure used a Nebulosa kernel-density plot, plot_eregulon_umaps() in grn/analysis/scenicplus_hybrid_all_template.qmd; it spreads
# the sparse signal -- MAFB is in ~10% of cells -- and gave little contrast.)
#
# Input: the hybrid-labeled Seurat object (common.R).
# Output: FIG4_OUT/panel_C_MAFB_umap.{pdf,png,svg}
#
#   Rscript figures/figure4/panel_C_mafb_umap.R
source(here::here("figures/figure4/common.R"))
suppressMessages({library(Seurat); library(qs2)})
GENE <- "MAFB"; SIZE <- 2.2
SEQ <- c("#dbe9f6", "#6baed6", "#2171b5", "#08306b")          # single-hue sequential (light -> dark blue)
fig_check_font()

obj <- qs_read(SEURAT_HYBRID_OBJ, nthreads = 8); DefaultAssay(obj) <- "SCT"
stopifnot(UMAP_REDUCTION %in% Reductions(obj), GENE %in% rownames(obj))
xy <- Embeddings(obj, UMAP_REDUCTION)[, 1:2]
d <- tibble(x = xy[, 1], y = xy[, 2], expr = FetchData(obj, vars = GENE, layer = "data")[, 1])
cap <- quantile(d$expr[d$expr > 0], 0.99)
cat(sprintf("%s: expressed in %.1f%% of %d cells; color scale capped at %.2f (99th percentile of expressing cells)\n", GENE, 100 * mean(d$expr > 0), nrow(d), cap))

p <- ggplot() +
  ggrastr::rasterise(geom_point(data = d %>% filter(expr == 0), aes(x, y), color = "#e3e2dc", size = 0.15, stroke = 0), dpi = 600) +
  ggrastr::rasterise(geom_point(data = d %>% filter(expr > 0) %>% arrange(expr), aes(x, y, color = pmin(expr, cap)), size = 0.25, stroke = 0), dpi = 600) +
  scale_color_gradientn(colors = SEQ, name = "log expr.", limits = c(0, cap), breaks = seq(0, floor(cap))) +
  coord_equal() +
  theme_void() + theme(legend.position = "bottom", legend.key.height = unit(0.07, "in"), legend.key.width = unit(0.35, "in"),
                       legend.title = element_text(size = FIG_PT_AXIS_TEXT, family = FIG_FONT), legend.text = element_text(size = FIG_PT_AXIS_TEXT, family = FIG_FONT),
                       plot.title = element_text(size = FIG_PT_AXIS_TITLE, family = FIG_FONT, hjust = 0.5), plot.background = element_blank(), panel.background = element_blank())
save_gg(p, "panel_C_MAFB_umap", width = SIZE, height = SIZE + 0.35)
