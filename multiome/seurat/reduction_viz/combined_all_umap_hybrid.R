library(Seurat)
library(tidyverse)
library(qs2)
library(seriation)
library(cowplot)
library(scCustomize)
library(scales)
library(here)
theme_set(theme_cowplot())

source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "scRNA.R"))
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "common_aesthetics.R"))

pnames = names(position_colors)
names(position_colors) = case_when(pnames == "nido" ~ "nr",
                                     pnames == "ncl" ~ "nc",
                                     TRUE ~ pnames)
# Directories -------------------------------------------------------------

## Parallels combined_all_umap.R, but on the hybrid snRNA-transferred labels from
## multiome/seurat/label_transfer_hybrid.qmd (`cluster_hybrid`) instead of the raw `cluster` identity
## that notebook curated by hand under the old naming scheme. qs2, not qs, because that notebook's
## output object is already qs2 -- this script reads/writes that format throughout rather than
## round-tripping through the older one.
naming_dir = here::here("multiome/seurat/label_transfer_hybrid")

data_fname = file.path(naming_dir, "obj_clustered_hybrid.qs2")
script_name = "combined_all_umap_hybrid"
out_dir = here::here("multiome/seurat/reduction_viz", script_name)
dir.create(out_dir, recursive = T)

data_out_obj_fname = file.path(out_dir, "obj_clustered.qs2")

# Load data ---------------------------------------------------------------

res_to_use = "cluster_hybrid"

dims_list = seq(30,50,10)
n.neighbors_list = c(30,40,50)
min.dist_list = c(0.1, 0.3, 0.5)

params = expand_grid(dims_list, n.neighbors_list, min.dist_list)

print(params)
set.seed(10)

redo = T
if (redo) {
  obj_filt = qs_read(data_fname, nthreads = 8)

  ## label_transfer_hybrid.qmd excludes Glut-Nido-3 (confirmed artefactual) before the transfer, so
  ## cluster_hybrid should already be NA-free -- asserted rather than assumed, so a future change to
  ## that notebook that reintroduces unlabeled cells fails loudly here instead of silently fitting the
  ## embedding on them.
  stopifnot("cluster_hybrid has NA cells -- check label_transfer_hybrid.qmd's exclusions" =
              !anyNA(obj_filt@meta.data[[res_to_use]]))

  DefaultAssay(obj_filt) = "SCT"
  obj_filt = RunPCA(obj_filt)

  # Integrate PCA embeddings
  obj_filt <- IntegrateLayers(object = obj_filt,
                              method = HarmonyIntegration,
                              orig.reduction = "pca",
                              new.reduction = 'pca_harmony', verbose = FALSE)

  for (i in 1:nrow(params)) {
    dims = params$dims_list[i]
    n.neighbors = params$n.neighbors_list[i]
    min.dist = params$min.dist_list[i]
    reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)

    print(reduction.name)
    obj_filt = RunUMAP(obj_filt,
                       reduction = "pca_harmony",
                           dims = 1:dims,
                           min.dist = min.dist,
                           n.neighbors = n.neighbors,
                           reduction.name=reduction.name
    )
  }

  qs_save(obj_filt, data_out_obj_fname, nthreads = 1)
} else {
  obj_filt = qs_read(data_out_obj_fname, nthreads = 8)
}

# Plot UMAP --------------------------------------------------------------------


reductions = Reductions(obj_filt)
reductions = grep("dims", reductions, value=T)
cats = c("position", res_to_use, "assignment")

for (reduction.name in reductions) {
  for ( ca in cats ) {
    ncat = length(unique(obj_filt@meta.data[,ca]))
    print(ncat)
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
    gg
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

  gg
  save_plot(file.path(out_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1 )
  save_plot(file.path(out_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=7, base_asp =1 )

  ca = res_to_use
  ncat = length(unique(obj_filt@meta.data[,ca]))
  gg = DimPlot_scCustom(obj_filt, reduction=reduction.name, group.by=ca, label=F, repel = T,
                        DiscretePalette_scCustomize(num_colors = ncat,
                                                    palette = "varibow")) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position="none"
    )

  gg
  save_plot(file.path(out_dir, sprintf("umap_%s_%s_no-label.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1 )
  save_plot(file.path(out_dir, sprintf("umap_%s_%s_no-label.png", reduction.name, ca)), gg, base_height=7, base_asp =1 )
}
