library(Seurat)
library(tidyverse)
library(qs)
library(seriation)
#library(proxy)
library(cowplot)
library(scCustomize)
library(scales)
theme_set(theme_cowplot())

source("../../../scripts/scRNA.R")
source("../../../scripts/common_aesthetics.R")

pnames = names(position_colors)
names(position_colors) = case_when(pnames == "nido" ~ "nr",
                                     pnames == "ncl" ~ "nc",
                                     TRUE ~ pnames)
#names(position_colors) = toupper(names(position_colors))
# Directories -------------------------------------------------------------

dir_root = "/ssd/brad/rstudio/multiome/motor-pathway/seurat/"
scrna_dir =  file.path(dir_root, "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/")

data_fname = file.path(scrna_dir, "obj_clustered.qs")
out_dir = file.path(scrna_dir, "reduction_viz")
script_name = "combined_all_umap"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive = T)

data_out_obj_fname = file.path(out_dir, "obj_clustered.qs")

# Load data ---------------------------------------------------------------

res_to_use = "cluster"

dims_list = seq(30,50,10)
n.neighbors_list = c(30,40,50)
min.dist_list = c(0.1, 0.3, 0.5)

params = expand_grid(dims_list, n.neighbors_list, min.dist_list)

print(params)
set.seed(10)

redo = T
if (redo) {
  obj_filt = qread(data_fname)
  DefaultAssay(obj_filt) = "SCT"
  obj_filt = RunPCA(obj_filt)
  
  # Integrate PCA embeddings
  obj_filt <- IntegrateLayers(object = obj_filt, 
                              method = HarmonyIntegration, 
                              orig.reduction = "pca", 
                              new.reduction = 'pca_harmony', verbose = FALSE)
  
  # obj_filt =  RunUMAP(
  #   object = obj_filt,
  #   reduction = "pca_harmony",
  #   reduction.name = "umap_rna_int", 
  #   assay = "SCT",
  #   dims = 1:rna_dims,
  #   verbose = TRUE
  # )
  
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
  
  qsave(obj_filt, data_out_obj_fname)
} else {
  obj_filt = qread(data_out_obj_fname)
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
