library(Seurat)
library(tidyverse)
library(qs)
library(seriation)
library(cowplot)
library(scCustomize)
library(scales)
library(Nebulosa)
theme_set(theme_cowplot())

source("../../scripts/scRNA.R")
source("../../scripts/common_aesthetics.R")
source("../../scripts/gene_lists.R")

pnames = names(position_colors)
names(position_colors) = case_when(pnames == "nido" ~ "nr",
                                     pnames == "ncl" ~ "nc",
                                     TRUE ~ pnames)
#names(position_colors) = toupper(names(position_colors))
# Directories -------------------------------------------------------------

dir_root = "/hdd/brad/rstudio/snRNA/snrna_cellranger/"
scrna_dir =  file.path(dir_root, "snrna_seurat_cellbender.0.01_preprocess")

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

redo = F
if (redo) {
  obj_int_filt = qread(data_fname)
  DefaultAssay(obj_int_filt) = "SCT"
  obj_int_filt = RunPCA(obj_int_filt)
  
  for (i in 1:nrow(params)) {
    dims = params$dims_list[i]
    n.neighbors = params$n.neighbors_list[i]
    min.dist = params$min.dist_list[i]
    reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)
    
    print(reduction.name)
    obj_int_filt = RunUMAP(obj_int_filt,
                           dims = 1:dims,
                           min.dist = min.dist,
                           n.neighbors = n.neighbors,
                           reduction.name=reduction.name
    )
  }
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

# Plot UMAP --------------------------------------------------------------------


reductions = Reductions(obj_int_filt)
reductions = grep("dims", reductions, value=T)
cats = c("position", res_to_use)

for (reduction.name in reductions) {
  for ( ca in cats ) {
    ncat = length(unique(obj_int_filt@meta.data[,ca]))
    gg = DimPlot_scCustom(obj_int_filt, 
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
  gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=F, repel = T ) +
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
  gg = DimPlot_scCustom(obj_int_filt, reduction=reduction.name, group.by=ca, label=F, repel = T, 
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


# Gene expression ---------------------------------------------------------

features = c(unlist(ct_markers), "SLC17A6", "GAD1", "SLC1A2", "CSF1R", "MBP", "PDGFRA", "NECTIN3", "SOX4", "FLI1",
             "SPEF2", "RGS5", "LUM", "NR4A2", "SIX2", "NDNF", "PVALB", "PCP4", "EMX2", "BACH2", "SOX2", "ISL1", "TBR1",
             "EOMES", "FOXP2", "MAFB", "LHX6", "PROX1", "LHX9", "NR2F2", "DACH1", "POU3F2",
             "SATB2", "SATB1", "FEZF2", "BCL11B", "EMX1", "ZBTB20")


plot_density_features <- function(obj, features, reduction, out_dir) {
  features_present <- intersect(features, rownames(obj))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  walk(features_present, function(feat) {
    gg <- plot_density(obj, features = feat, reduction = reduction, size = 0.1) +
      coord_equal() + umap_theme()
    print(gg)
    save_plot(file.path(out_dir, sprintf("umap_density_%s.pdf", feat)), gg,
              base_height = 5, base_width = 5)
  })
}

reduction = "dims40nn30mindist0.3"
density_dir = file.path(out_dir, "density")
plot_density_features(obj_int_filt, features, reduction = reduction, out_dir = density_dir)


