library(Seurat)
library(tidyverse)
library(qs)
library(seriation)
library(cowplot)
library(scCustomize)
library(scales)
theme_set(theme_cowplot())

source("../../scripts/scRNA.R")
source("../../scripts/common_aesthetics.R")

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
script_name = "combined_all_umap_glut"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive = T)

data_out_obj_fname = file.path(out_dir, "obj_clustered.qs")

# Load data ---------------------------------------------------------------

res_to_use = "cluster"

dims_list = seq(10,30,10)
n.neighbors_list = c(30,40,50)
min.dist_list = c(0.3, 0.4, 0.5)

params = expand_grid(dims_list, n.neighbors_list, min.dist_list)

print(params)
set.seed(10)

redo = F
if (redo) {
  obj_int_filt = qread(data_fname)
  obj_int_filt$region = case_when(obj_int_filt$position %in% c("arco", "ra") ~ "arco",
                                  obj_int_filt$position %in% c("nc", "hvc", "nr", "lman") ~ "nido")
  DefaultAssay(obj_int_filt) = "SCT"
  
  cells = Cells(obj_int_filt)[grepl("Glut", obj_int_filt$cluster)]
  obj_int_filt = subset(obj_int_filt, cells=cells)
  
  cells_arco = Cells(obj_int_filt)[grepl("RA|Arco", obj_int_filt$cluster)]
  cells_nido = Cells(obj_int_filt)[grepl("HVC|NC|Nido|LMAN|NR|Pre", obj_int_filt$cluster)]
  objs = list(arco = subset(obj_int_filt, cells=cells_arco),
              nido = subset(obj_int_filt, cells=cells_nido))
  
  objs = map(objs, function(obj_int_filt) {
    obj_int_filt = obj_int_filt |>
      SCTransform() |>
      RunPCA()
    
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
    obj_int_filt
  })
  
  qsave(objs, data_out_obj_fname)
} else {
  objs = qread(data_out_obj_fname)
}

# Plot UMAP --------------------------------------------------------------------

iwalk(objs, function(obj_int_filt, region_cur) {
  reductions = Reductions(obj_int_filt)
  reductions = grep("dims", reductions, value=T)
  cats = c("position", res_to_use)
  
  for (reduction.name in reductions) {
    for ( ca in cats ) {
      ncat = length(unique(obj_int_filt@meta.data[,ca]))
      print(ncat)
      
      # if (ncat < 10)  {
      #   pal = "stepped"
      # } else {
        pal = "varibow"
      #}
      gg = DimPlot_scCustom(obj_int_filt, 
                            reduction=reduction.name, group.by=ca, label=T, repel = T, pt.size = 2,
                            raster = T, raster.dpi = c(1024,1024),
                            DiscretePalette_scCustomize(num_colors = ncat,
                                                        palette = pal)) +
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position="none"
        ) +
        labs(x="", y="UMAP2") 
      gg
      save_plot(file.path(out_dir, sprintf("umap_%s_%s_%s.pdf", region_cur, ca, reduction.name)), gg, base_height=7, base_asp =1)
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
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_%s.pdf", region_cur,  reduction.name, ca)), gg, base_height=7, base_asp =1 )
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_%s.png", region_cur, reduction.name, ca)), gg, base_height=7, base_asp =1 )
    
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
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_%s_no-label.pdf", region_cur, reduction.name, ca)), gg, base_height=7, base_asp =1 )
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_%s_no-label.png", region_cur, reduction.name, ca)), gg, base_height=7, base_asp =1 )
  }
})