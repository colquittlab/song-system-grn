library(Seurat)
library(Signac)
library(tidyverse)
library(qs2)
library(cowplot)
library(scCustomize)
library(scales)
library(here)
theme_set(theme_cowplot())

source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "scRNA.R"))
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "common_aesthetics.R"))
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))

pnames = names(position_colors)
names(position_colors) = case_when(pnames == "nido" ~ "nr",
                                     pnames == "ncl" ~ "nc",
                                     TRUE ~ pnames)
# Directories -------------------------------------------------------------

## Hybrid-label counterpart of the pre-hybrid combined_glut_umap_dotplot_tfs.R (motor-pathway
## reduction_viz), whose glut-only object multiome/archr/pos_regulators.qmd plotted its UMAPs on.
## That script kept every Glut cluster except the precursors ("Pre"); under the hybrid names that is
## hybrid_groups$glut minus hybrid_groups$precursor. The reductions (PCA, Harmony, LSI, WNN, UMAP)
## are all recomputed on that subset; the SCT normalization is inherited from the full object, as it
## was in the original, so expression values match the other hybrid figures.
data_fname = HYBRID_OBJ_QS2
script_name = "combined_glut_umap_hybrid"
out_dir = here::here("multiome/seurat/reduction_viz", script_name)
dir.create(out_dir, recursive = T)

data_out_obj_fname = file.path(out_dir, "obj_clustered.qs2")

# Load data ---------------------------------------------------------------

res_to_use = "cluster_hybrid"

## The original swept 20/30/40 dims for each modality; pos_regulators.qmd used rna30/atac30 only.
rna_dims_list = c(30)
atac_dims_list = c(30)
n.neighbors_list = c(30)
min.dist_list = c(0.3)

set.seed(10)

redo = T
if (redo) {
  obj = qs_read(data_fname, nthreads = 8)
  glut_clusters = setdiff(hybrid_groups$glut, hybrid_groups$precursor)
  cells = Cells(obj)[obj$cluster_hybrid %in% glut_clusters]
  obj_filt = subset(obj, cells=cells)
  rm(obj); gc()

  ## Only RNA (for the per-run split), ATAC and SCT are used downstream; drop the rest, and the
  ## full-object reductions so nothing computed on all cells can leak into the glut embedding.
  DefaultAssay(obj_filt) = "SCT"
  obj_filt = DietSeurat(obj_filt, assays = c("RNA", "ATAC", "SCT"), layers = NULL,
                        dimreducs = character(), graphs = character())

  obj_filt = process_multimodal_pipeline(obj_filt,
                                         batch_var="run",
                                         cluster_prefix = "cluster_glut",
                                         rna_dims_list = rna_dims_list,
                                         atac_dims_list = atac_dims_list,
                                         min.dist_list = min.dist_list,
                                         n.neighbors_list = n.neighbors_list,
                                         run_harmony = T,
                                         run_sct = T,
                                         umap_prefix = "umap",
                                         run_parallel = F,
                                         run_clusters = F
  )
  qs_save(obj_filt, data_out_obj_fname, nthreads = 8)
} else {
  obj_filt = qs_read(data_out_obj_fname, nthreads = 8)
}

# Plot UMAP --------------------------------------------------------------------


reductions = Reductions(obj_filt)
reductions = grep("umap", reductions, value=T)
cats = c("position", res_to_use, "assignment")

for (reduction.name in reductions) {
  for ( ca in cats ) {
    ncat = length(unique(obj_filt@meta.data[,ca]))
    gg = DimPlot_scCustom(obj_filt,
                          reduction=reduction.name, group.by=ca, label=T, repel = T, pt.size = 2,
                          raster = T, raster.dpi = c(1024,1024),
                          DiscretePalette_scCustomize(num_colors = ncat,
                                                      palette = "varibow")) +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position="none"
      ) +
      labs(x="", y="UMAP2")
    save_plot(file.path(out_dir, sprintf("umap_%s_%s.pdf", ca, reduction.name)), gg, base_height=7, base_asp =1)
  }

  ca = "position"
  gg = DimPlot(obj_filt, reduction=reduction.name, group.by=ca, label=F, repel = T ) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position="none"
    ) +
    scale_color_manual(values=position_colors)
  save_plot(file.path(out_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1 )
  save_plot(file.path(out_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=7, base_asp =1 )
}
