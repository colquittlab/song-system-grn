#!/usr/bin/env Rscript
# Figure 4, panel C -- MAFB expression on the multiome UMAP.
#
# Kernel-density estimate of MAFB expression (Nebulosa::plot_density, SCT log-normalized data) on the hybrid-labeled object's UMAP (reduction
# dims25nn30mindist0.3: Harmony-integrated PCA, 25 dimensions, 30 neighbours, min.dist 0.3; multiome/seurat/reduction_viz/combined_all_umap_hybrid.R).
# Same call as plot_eregulon_umaps() in grn/analysis/scenicplus_hybrid_all_template.qmd (rendered by grn/analysis/all.qmd), which wrote
#   <scenicplus output>/config37/umaps/MAFB_extended_+:+/MAFB_density.pdf
#
# Input: the hybrid-labeled Seurat object (common.R).
# Output: FIG4_OUT/panel_C_MAFB_density.pdf (and the non-smoothed feature plot, panel_C_MAFB.pdf)
#
#   Rscript figures/figure4/panel_C_mafb_umap_density.R
source(here::here("figures/figure4/common.R"))
suppressMessages({library(Seurat); library(scCustomize); library(qs2); library(cowplot)})
GENE <- "MAFB"
obj <- qs_read(SEURAT_HYBRID_OBJ, nthreads = 8)
DefaultAssay(obj) <- "SCT"
stopifnot(UMAP_REDUCTION %in% Reductions(obj), GENE %in% rownames(obj))
gg <- Nebulosa::plot_density(obj, features = GENE, size = .1, reduction = UMAP_REDUCTION) + coord_equal()
save_plot(file.path(FIG4_OUT, paste0("panel_C_", GENE, "_density.pdf")), gg, base_width = 7, base_asp = 1)
gg <- FeaturePlot_scCustom(obj, features = GENE, figure_plot = TRUE, aspect_ratio = 1, reduction = UMAP_REDUCTION, raster = TRUE, raster.dpi = c(1024, 1024),
                           pt.size = 2, alpha_exp = 1, alpha_na_exp = .1)
save_plot(file.path(FIG4_OUT, paste0("panel_C_", GENE, ".pdf")), gg, base_width = 7, base_asp = 1)
message("saved ", file.path(FIG4_OUT, paste0("panel_C_", GENE, "_density.pdf")))
