library(Seurat)
library(tidyverse)
library(qs2)
library(seriation)
library(cowplot)
library(scCustomize)
library(scales)
library(Nebulosa)
library(here)
theme_set(theme_cowplot())

source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "scRNA.R"))
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "common_aesthetics.R"))
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "gene_lists.R"))

pnames = names(position_colors)
names(position_colors) = case_when(pnames == "nido" ~ "nr",
                                     pnames == "ncl" ~ "nc",
                                     TRUE ~ pnames)
# Directories -------------------------------------------------------------

## Parallels combined_all_umap.R, but on the current (hybrid division-based) labels from
## snrna/naming/hybrid_division_naming.qmd instead of the raw `cluster` identity, and with the
## clusters that notebook excludes (Glut-Arco-4, Glut-Nido-3 -- confirmed artifacts) dropped from
## the object entirely before the embedding is fit, not merely left unlabeled on top of it. qs2, not
## qs, because that notebook's output object is already qs2 -- this script reads/writes that format
## throughout rather than round-tripping through the older one.
naming_dir = here::here("snrna/naming/hybrid_division_naming")

data_fname = file.path(naming_dir, "obj_hybrid_labels.qs2")
script_name = "combined_all_umap_hybrid"
out_dir = here::here("snrna/reduction_viz", script_name)
dir.create(out_dir, recursive = T)

data_out_obj_fname = file.path(out_dir, "obj_clustered.qs2")

# Load data ---------------------------------------------------------------

res_to_use = "celltype_hybrid"

dims_list = seq(30,50,10)
n.neighbors_list = c(30,40,50)
min.dist_list = c(0.1, 0.3, 0.5)

params = expand_grid(dims_list, n.neighbors_list, min.dist_list)

print(params)
set.seed(10)

redo = T
if (redo) {
  obj_int_filt = qs_read(data_fname, nthreads = 8)

  ## Cluster omissions: celltype_hybrid is NA for exactly the clusters that notebook excluded
  ## (Glut-Arco-4, Glut-Nido-3) -- dropped here so the embedding itself is never fit on them, not
  ## just left out of the color legend.
  obj_int_filt = subset(obj_int_filt, cells = colnames(obj_int_filt)[!is.na(obj_int_filt@meta.data[[res_to_use]])])

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

  qs_save(obj_int_filt, data_out_obj_fname, nthreads = 8)
} else {
  obj_int_filt = qs_read(data_out_obj_fname, nthreads = 8)
}

# Individual identity (souporcell) ----------------------------------------

## Bird identity is not in the object as loaded. `assignment` is souporcell's *per-library* cluster
## index (0/1/2), recycled independently in each of the six libraries, so it names three arbitrary
## labels six times over rather than three birds; snrna/clustering/snrna_souporcell_clustering.qmd
## replaces it by matching clusters across libraries on genotype. That notebook's per-cell table is
## read here rather than re-derived, and joined on cell name -- which is exact: every cell of this
## object is in that table (asserted below), both being built from the same libraries in the same
## order, so the `-<sample_index>` suffix means the same thing on both sides.
##
## Two properties of that result constrain every plot below.
##
## 1. There are six birds in two batches that share no individual: bird1-3 in ra/arco/hvc/nc,
##    bird4-6 in lman/nr. A single "colour by individual" panel over all cells is therefore six
##    colours, not three, and two colours from different batches say nothing about each other --
##    hence per-batch panels alongside the combined one.
## 2. `position` is not the library. The preprocessing notebook re-derives position for the arco
##    library from cluster identity (Glut-RA -> "ra"), so ~500 cells sequenced in ra_run1 carry
##    position "arco" and 22 arco cells carry "ra". Batch is a property of the *library*, so it is
##    taken from souporcell's `sample_id` (which agrees with `sample_name` cell-for-cell), never
##    from `position`.
soup_fname = here::here("snrna/clustering/snrna_souporcell_clustering", "souporcell_assignments.csv")

add_individual_metadata = function(obj) {
  soup = read_csv(soup_fname, show_col_types = FALSE)

  stopifnot(
    "souporcell assignments do not cover every cell -- check the cell-name suffix convention" =
      all(colnames(obj) %in% soup$cell)
  )

  md = tibble(cell = colnames(obj)) %>%
    left_join(soup %>% select(cell, sample_id, individual, batch, cluster),
              by = "cell") %>%
    mutate(
      individual = factor(individual, levels = sprintf("bird%s", 1:6)),
      soup_batch = factor(sprintf("batch%s", batch)),
      soup_library = sample_id,
      soup_cluster = cluster
    ) %>%
    select(cell, individual, soup_batch, soup_library, soup_cluster) %>%
    column_to_rownames("cell")

  ## Every cell here is already a souporcell singlet (the preprocessing notebook dropped the
  ## doublets), so `individual` is never NA in practice -- assert rather than assume, because a
  ## silent NA would plot as a seventh "bird".
  stopifnot("unassigned individuals after the join" = !any(is.na(md$individual)))
  AddMetaData(obj, md)
}

if (!"individual" %in% colnames(obj_int_filt@meta.data)) {
  obj_int_filt = add_individual_metadata(obj_int_filt)
  ## Write the identity back into the saved object so downstream scripts -- notably
  ## snrna/variation/cross_individual_variation.R -- read it rather than repeating the join.
  qs_save(obj_int_filt, data_out_obj_fname, nthreads = 8)
}

print(table(obj_int_filt$individual, obj_int_filt$soup_library))

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


# UMAP by individual ------------------------------------------------------

## One reduction only -- dims40nn30mindist0.3, the same one the density panels below use.
##
## Plotted from the embedding by hand rather than through DimPlot(split.by=) for two reasons. The
## draw order matters: whichever bird is drawn last covers the others, and with birds occupying the
## same cell types that is the difference between "bird3 is everywhere" and the truth, so the points
## are shuffled before drawing. And the split panels need the *other* birds' cells behind them in
## grey -- without that background each panel is a different-shaped point cloud and nothing in it is
## comparable across panels.
reduction_ind = "dims40nn30mindist0.3"

individual_levels = levels(obj_int_filt$individual)
individual_colors = setNames(
  DiscretePalette_scCustomize(num_colors = length(individual_levels), palette = "varibow"),
  individual_levels
)

## Which libraries a batch was sequenced from -- from the data, not restated, so a re-run of the
## souporcell notebook that changes the grouping changes these panels too.
batch_libraries = obj_int_filt@meta.data %>%
  count(soup_batch, soup_library) %>%
  group_by(soup_batch) %>%
  summarise(libraries = paste(sort(unique(soup_library)), collapse = ", "), .groups = "drop")
print(batch_libraries)

umap_df = function(obj, reduction, vars) {
  emb = Embeddings(obj, reduction)[, 1:2]
  colnames(emb) = c("UMAP_1", "UMAP_2")
  bind_cols(as_tibble(emb, rownames = "cell"),
            obj@meta.data[, vars, drop = FALSE] %>% as_tibble())
}

ind_df = umap_df(obj_int_filt, reduction_ind,
                 c("individual", "soup_batch", "soup_library", "position", res_to_use))

## theme_cowplot leaves plot.background transparent, which is invisible against a white viewer and
## black against a dark one -- fine for the solid-colour panels above, not for these, where the
## point of the figure is grey cells against the background. Paint it white explicitly.
umap_theme_white = function() {
  umap_theme() + theme(plot.background = element_rect(fill = "white", colour = NA))
}

plot_by_individual = function(df, pt.size = 0.3) {
  ggplot(df[sample(nrow(df)), ], aes(UMAP_1, UMAP_2, colour = individual)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_colour_manual(values = individual_colors, drop = FALSE, name = NULL) +
    coord_equal() +
    umap_theme_white() +
    guides(colour = guide_legend(override.aes = list(size = 3)))
}

split_by_individual = function(df, ncol = 3, pt.size = 0.3) {
  ggplot(df, aes(UMAP_1, UMAP_2)) +
    ## Background: every cell in the panel set, with the facet column dropped so it repeats behind
    ## each panel.
    geom_point(data = df %>% select(-individual), colour = "grey85", size = pt.size, stroke = 0) +
    geom_point(aes(colour = individual), size = pt.size, stroke = 0) +
    scale_colour_manual(values = individual_colors, drop = FALSE, guide = "none") +
    facet_wrap(~ individual, ncol = ncol) +
    coord_equal() +
    umap_theme_white()
}

gg = plot_by_individual(ind_df)
save_plot(file.path(out_dir, sprintf("umap_%s_individual.pdf", reduction_ind)), gg, base_height = 7, base_asp = 1.15)
save_plot(file.path(out_dir, sprintf("umap_%s_individual.png", reduction_ind)), gg, base_height = 7, base_asp = 1.15)

gg = split_by_individual(ind_df, ncol = 3)
save_plot(file.path(out_dir, sprintf("umap_%s_individual_split.pdf", reduction_ind)), gg, base_height = 7, base_asp = 1.5)
save_plot(file.path(out_dir, sprintf("umap_%s_individual_split.png", reduction_ind)), gg, base_height = 7, base_asp = 1.5)

## Per batch. Birds are only comparable within a batch, and the two batches cover different regions
## (bird1-3 ra/arco/hvc/nc, bird4-6 lman/nr), so the combined split panels above put a
## four-region cloud next to a two-region one. These restrict the background to the batch's own
## cells, which is the comparison that is actually meaningful.
for (b in levels(droplevels(ind_df$soup_batch))) {
  df_b = ind_df %>% filter(soup_batch == b) %>% mutate(individual = droplevels(individual))

  gg = plot_by_individual(df_b)
  save_plot(file.path(out_dir, sprintf("umap_%s_individual_%s.pdf", reduction_ind, b)), gg, base_height = 7, base_asp = 1.15)
  save_plot(file.path(out_dir, sprintf("umap_%s_individual_%s.png", reduction_ind, b)), gg, base_height = 7, base_asp = 1.15)

  gg = split_by_individual(df_b, ncol = 3)
  save_plot(file.path(out_dir, sprintf("umap_%s_individual_split_%s.pdf", reduction_ind, b)), gg, base_height = 4, base_asp = 2.6)
  save_plot(file.path(out_dir, sprintf("umap_%s_individual_split_%s.png", reduction_ind, b)), gg, base_height = 4, base_asp = 2.6)
}

## Library alongside individual, because in this design the two are confounded in opposite
## directions: a bird spans several libraries, a library holds three birds, and a bird-shaped patch
## of the embedding that is also library-shaped is a library effect. Faceting bird x library makes
## which of the two is driving a patch visible in one figure.
gg = ggplot(ind_df, aes(UMAP_1, UMAP_2)) +
  geom_point(data = ind_df %>% select(-individual, -soup_library), colour = "grey88", size = 0.2, stroke = 0) +
  geom_point(aes(colour = individual), size = 0.2, stroke = 0) +
  scale_colour_manual(values = individual_colors, drop = FALSE, guide = "none") +
  facet_grid(soup_library ~ individual) +
  coord_equal() +
  umap_theme_white()
## png only: 36 panels of the full embedding is a 60 MB vector file and nothing in it is read at
## point level.
save_plot(file.path(out_dir, sprintf("umap_%s_individual_by_library.png", reduction_ind)), gg, base_height = 12, base_asp = 0.85)

## Cells per bird per cell type -- the count behind every panel above, and the sample sizes the
## cross-individual expression analysis is limited by.
ind_df %>%
  count(soup_batch, soup_library, individual, .data[[res_to_use]]) %>%
  write_csv(file.path(out_dir, "cells_per_individual_celltype.csv"))

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

