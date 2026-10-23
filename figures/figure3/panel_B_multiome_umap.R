#!/usr/bin/env Rscript
# Figure 3, panel B -- UMAP of the hybrid-labeled multiome object, all cell types, colored by cell type (cluster_hybrid).
#
# Same drawing as the "_no-label" plot of multiome/seurat/reduction_viz/combined_all_umap_hybrid.R (the version that made the figure files): the saved
# reduction dims25nn30mindist0.3 (Harmony-integrated PCA, 25 dimensions, 30 neighbors, min.dist 0.3), drawn with scCustomize::DimPlot_scCustom
# (default point size and shuffled drawing order, seed 1), the "varibow" palette with one color per cell type, axes, legend and labels removed,
# 7 x 7 in (cowplot::save_plot). The UMAP is read from the object, not recomputed.
# Everything else in the panel was drawn by hand in Inkscape (cell-type labels, gray shaded regions, UMAP axis arrows, "n = 16,588"); see README.md.
#
# Input: the hybrid-labeled Seurat object (common.R).
# Output: FIG3_OUT/panel_B_umap_dims25nn30mindist0.3_cluster_hybrid_no-label.{pdf,png} (+ panel_B_cluster_colors.csv: cell type, color, n cells)
#
#   Rscript figures/figure3/panel_B_multiome_umap.R
source(here::here("figures/figure3/common.R"))
suppressMessages({library(Seurat); library(qs2); library(scCustomize); library(cowplot)})

UMAP_REDUCTION <- "dims25nn30mindist0.3"; GROUP <- "cluster_hybrid"
NAME <- sprintf("panel_B_umap_%s_%s_no-label", UMAP_REDUCTION, GROUP)
# as in the source script: the plot's theme elements not set by DimPlot are filled from the session theme, so cowplot's (7 pt plot margin) sets the
# panel size; with ggplot2's default theme the panel is 0.6 % larger
theme_set(theme_cowplot())

obj <- qs_read(SEURAT_HYBRID_OBJ, nthreads = 8)
stopifnot(UMAP_REDUCTION %in% Reductions(obj), !anyNA(obj@meta.data[[GROUP]]))
ncat <- length(unique(obj@meta.data[[GROUP]]))
pal <- DiscretePalette_scCustomize(num_colors = ncat, palette = "varibow")

gg <- DimPlot_scCustom(obj, reduction = UMAP_REDUCTION, group.by = GROUP, label = FALSE, repel = TRUE, colors_use = pal) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "none")

lv <- levels(factor(obj@meta.data[[GROUP]]))
cols <- tibble(cluster = lv, color = unname(pal[seq_along(lv)]), n_cells = as.integer(table(factor(obj@meta.data[[GROUP]], levels = lv))))
cat(sprintf("%s nuclei, %d cell types\n", format(ncol(obj), big.mark = ","), ncat))
print(as.data.frame(cols))
write_csv(cols, file.path(FIG3_OUT, "panel_B_cluster_colors.csv"))

save_plot(file.path(FIG3_OUT, paste0(NAME, ".pdf")), gg, base_height = 7, base_asp = 1)
save_plot(file.path(FIG3_OUT, paste0(NAME, ".png")), gg, base_height = 7, base_asp = 1)
message("saved ", file.path(FIG3_OUT, NAME))
