#!/usr/bin/env Rscript
# Figure 1, panels D and E -- UMAP of the snRNA-seq data coloured by cell type (D) and by dissected region (E).
#
# PCA on the SCT assay of the labelled cells (Seurat::RunPCA defaults), UMAP on PCs 1-40 with 30 neighbours and min.dist 0.3 (Seurat::RunUMAP,
# uwot, seed 42). Same calls as snrna/reduction_viz/combined_all_umap_hybrid.R, which fits a grid of UMAPs and wrote, for this one
# (dims40nn30mindist0.3),
#   snrna/reduction_viz/combined_all_umap_hybrid/umap_dims40nn30mindist0.3_celltype_hybrid_no-label.png   (D)
#   snrna/reduction_viz/combined_all_umap_hybrid/umap_dims40nn30mindist0.3_position.png                   (E)
# Refitting from the labelled object reproduces that script's saved embedding exactly. Cell-type and region names in the figure were added by hand.
#
# Input: the snRNA-seq object with the hybrid labels (common.R).
# Output: FIG1_OUT/panel_D_umap_celltype.{pdf,png}, FIG1_OUT/panel_E_umap_position.{pdf,png}, FIG1_OUT/panel_DE_umap_embedding.csv
#
#   Rscript figures/figure1/panel_DE_snrna_umap.R
source(here::here("figures/figure1/common.R"))
suppressMessages({library(Seurat); library(scCustomize); library(cowplot)})
theme_set(theme_cowplot())
REDUCTION <- "dims40nn30mindist0.3"

obj <- read_snrna()
stopifnot("unexpected dissection positions" = setequal(unique(obj$position), POSITIONS))
set.seed(10)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:40, n.neighbors = 30, min.dist = 0.3, reduction.name = REDUCTION, verbose = FALSE)
as_tibble(Embeddings(obj, REDUCTION), rownames = "cell") %>%
  bind_cols(obj@meta.data[, c(CELLTYPE_COL, "position")]) %>%
  write_csv(file.path(FIG1_OUT, "panel_DE_umap_embedding.csv"))

umap_bare <- theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
                   legend.position = "none")
save_both <- function(gg, name) {
  for (ext in c("pdf", "png")) save_plot(file.path(FIG1_OUT, paste0(name, ".", ext)), gg, base_height = 7, base_asp = 1)
  message("saved ", file.path(FIG1_OUT, name))
}

# D: cell types, scCustomize "varibow" palette over the 47 types (in alphabetical order of the labels)
n_types <- length(unique(obj@meta.data[[CELLTYPE_COL]]))
gg <- DimPlot_scCustom(obj, reduction = REDUCTION, group.by = CELLTYPE_COL, label = FALSE, repel = TRUE,
                       colors_use = DiscretePalette_scCustomize(num_colors = n_types, palette = "varibow")) + umap_bare
save_both(gg, "panel_D_umap_celltype")

# E: dissected region
gg <- DimPlot(obj, reduction = REDUCTION, group.by = "position", label = FALSE, repel = TRUE) + umap_bare + scale_color_manual(values = POSITION_COLORS)
save_both(gg, "panel_E_umap_position")
