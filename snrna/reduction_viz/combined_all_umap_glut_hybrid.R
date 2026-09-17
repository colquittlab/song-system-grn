library(Seurat)
library(tidyverse)
library(qs2)
library(seriation)
library(cowplot)
library(scCustomize)
library(scales)
library(here)
theme_set(theme_cowplot())
options(future.globals.maxSize = Inf)

source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "scRNA.R"))
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "common_aesthetics.R"))

pnames = names(position_colors)
names(position_colors) = case_when(pnames == "nido" ~ "nr",
                                     pnames == "ncl" ~ "nc",
                                     TRUE ~ pnames)

# Directories -------------------------------------------------------------

## Parallels combined_all_umap_glut.R, but on the current (hybrid division-based) labels from
## snrna/naming/hybrid_division_naming.qmd instead of the raw `cluster` identity. As in
## combined_all_umap_hybrid.R, the clusters that notebook excludes (Glut-Arco-4, Glut-Nido-3 --
## confirmed artifacts) are dropped from the object before the embedding is fit, not merely left
## unlabeled on top of it, and the format is qs2 throughout because that notebook's output is qs2.
naming_dir = here::here("snrna/naming/hybrid_division_naming")

data_fname = file.path(naming_dir, "obj_hybrid_labels.qs2")
script_name = "combined_all_umap_glut_hybrid"
out_dir = here::here("snrna/reduction_viz", script_name)
dir.create(out_dir, recursive = T)

data_out_obj_fname = file.path(out_dir, "obj_clustered.qs2")

# Load data ---------------------------------------------------------------

res_to_use = "celltype_hybrid"

dims_list = seq(10,30,10)
n.neighbors_list = c(30,40,50)
min.dist_list = c(0.3, 0.4, 0.5)

params = expand_grid(dims_list, n.neighbors_list, min.dist_list)

print(params)
set.seed(10)

redo = T
if (redo) {
  obj_int_filt = qs_read(data_fname, nthreads = 8)
  obj_int_filt$region = case_when(obj_int_filt$position %in% c("arco", "ra") ~ "arco",
                                  obj_int_filt$position %in% c("nc", "hvc", "nr", "lman") ~ "nido")
  DefaultAssay(obj_int_filt) = "SCT"

  ## Cluster omissions: celltype_hybrid is NA for exactly the clusters the naming notebook excluded
  ## (Glut-Arco-4, Glut-Nido-3) -- dropped here so neither embedding is ever fit on them, rather
  ## than left out of the colour legend only.
  obj_int_filt = subset(obj_int_filt, cells = colnames(obj_int_filt)[!is.na(obj_int_filt@meta.data[[res_to_use]])])

  ## The excitatory set. Under the hybrid scheme every excitatory label still carries the `Glut-`
  ## prefix -- the divisions (Glut-CACNA1H-*, Glut-DACH2-*, Glut-SATB2-*), the song-region clusters
  ## (…-HVCra, …-HVCra-Int, …-HVCx, …-LMANsh, …-LMANco, …-RA), the precursor series
  ## (Glut-NSC/NB/Im) and Glut-GABA -- so this selects the cells the `cluster`-based script did,
  ## minus the two excluded clusters. GABA-Im (ex GABA-Pre) is correctly not caught, as before.
  cells = Cells(obj_int_filt)[grepl("^Glut", obj_int_filt@meta.data[[res_to_use]])]
  obj_int_filt = subset(obj_int_filt, cells=cells)

  ## Arcopallial vs nidopallial split. In the old labels this was a grepl over region-named clusters
  ## ("RA|Arco" against "HVC|NC|Nido|LMAN|NR|Pre"); the hybrid scheme names the arcopallial division
  ## CACNA1H outright, so arco is that division and nido is the rest of the excitatory set (DACH2,
  ## SATB2, the precursor series, Glut-GABA). Taking nido as the complement rather than as a second
  ## pattern means a label added to the scheme later cannot silently fall out of both halves.
  is_arco = grepl("CACNA1H", obj_int_filt@meta.data[[res_to_use]])
  cells_arco = Cells(obj_int_filt)[is_arco]
  cells_nido = Cells(obj_int_filt)[!is_arco]
  stopifnot("arco/nido split does not partition the excitatory cells" =
              length(cells_arco) + length(cells_nido) == ncol(obj_int_filt))
  print(table(obj_int_filt@meta.data[[res_to_use]], if_else(is_arco, "arco", "nido")))

  objs = list(arco = subset(obj_int_filt, cells=cells_arco),
              nido = subset(obj_int_filt, cells=cells_nido))

  objs = map(objs, function(obj_int_filt) {
    ## celltype_hybrid is a factor over all 47 labels of the full object; dropping the levels that
    ## are not in this half keeps the palette from being sized to labels no panel can show.
    obj_int_filt@meta.data[[res_to_use]] = droplevels(factor(obj_int_filt@meta.data[[res_to_use]]))

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

  qs_save(objs, data_out_obj_fname, nthreads = 8)
} else {
  objs = qs_read(data_out_obj_fname, nthreads = 8)
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

      pal = "varibow"
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
    ncat = length(unique(obj_int_filt@meta.data[,ca]))
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
