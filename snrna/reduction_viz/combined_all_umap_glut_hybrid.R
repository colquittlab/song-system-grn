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

  ## Further omissions, specific to this script rather than to the naming notebook: Glut-SATB2-1 is
  ## the lone member of its division, Glut-GABA is the mixed-identity cluster (ex Glut-Nido-2), and
  ## Glut-NSC is the neural stem cell head of the precursor series. All three sit far enough off the
  ## DACH2 body to set the scale of the nidopallial embedding while saying nothing about it, so they
  ## are dropped before the UMAP is fit, not just hidden after. All three are nidopallial, so the
  ## arcopallial half is unaffected. They remain in obj_hybrid_labels.qs2 and in every other
  ## consumer of it.
  umap_excluded_labels = c("Glut-SATB2-1", "Glut-GABA", "Glut-NSC")
  stopifnot("a umap_excluded_labels entry is not a celltype_hybrid label" =
              all(umap_excluded_labels %in% unique(obj_int_filt@meta.data[[res_to_use]])))
  obj_int_filt = subset(obj_int_filt,
                        cells = Cells(obj_int_filt)[!obj_int_filt@meta.data[[res_to_use]] %in% umap_excluded_labels])

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

# Individual identity (souporcell) ----------------------------------------

## Bird identity is not in the object as loaded. `assignment` is souporcell's *per-library* cluster
## index (0/1/2), recycled independently in each of the six libraries, so it names three arbitrary
## labels six times over rather than three birds; snrna/clustering/snrna_souporcell_clustering.qmd
## replaces it by matching clusters across libraries on genotype. That notebook's per-cell table is
## read here rather than re-derived, and joined on cell name, which is exact (asserted below).
##
## Two properties of that result constrain the panels: there are six birds in two batches that share
## no individual (bird1-3 in ra/arco/hvc/nc, bird4-6 in lman/nr), so two colours from different
## batches say nothing about each other; and batch is a property of the *library*, so it is taken
## from souporcell's `sample_id`, never from `position` -- the preprocessing notebook re-derives
## position for the arco library from cluster identity, so the two disagree for some hundreds of
## cells.
soup_fname = here::here("snrna/clustering/snrna_souporcell_clustering", "souporcell_assignments.csv")

add_individual_metadata = function(obj) {
  soup = read_csv(soup_fname, show_col_types = FALSE)

  stopifnot(
    "souporcell assignments do not cover every cell -- check the cell-name suffix convention" =
      all(colnames(obj) %in% soup$cell)
  )

  md = tibble(cell = colnames(obj)) %>%
    left_join(soup %>% select(cell, sample_id, individual, batch), by = "cell") %>%
    mutate(individual = factor(individual, levels = sprintf("bird%s", 1:6)),
           soup_batch = factor(sprintf("batch%s", batch)),
           soup_library = sample_id) %>%
    select(cell, individual, soup_batch, soup_library) %>%
    column_to_rownames("cell")

  ## Every cell here is already a souporcell singlet (the preprocessing notebook dropped the
  ## doublets), so `individual` is never NA in practice -- assert rather than assume, because a
  ## silent NA would plot as a seventh "bird".
  stopifnot("unassigned individuals after the join" = !any(is.na(md$individual)))
  AddMetaData(obj, md)
}

if (!all(map_lgl(objs, ~ "individual" %in% colnames(.x@meta.data)))) {
  objs = map(objs, add_individual_metadata)
  ## Write the identity back into the saved object so anything reading this embedding -- notably
  ## snrna/deg/song-surround_deg_glut_hybrid.qmd -- gets it without repeating the join.
  qs_save(objs, data_out_obj_fname, nthreads = 8)
}

iwalk(objs, function(obj_cur, region_cur) {
  cat("\n==", region_cur, "-- individual x library\n")
  print(table(obj_cur$individual, obj_cur$soup_library))
})

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

# UMAP by individual ------------------------------------------------------

## The same embeddings, coloured by bird instead of by label -- the panel that says whether a
## cluster that separates cleanly is separating by cell type or by donor. Written for every
## reduction, so each cell-type panel above has a bird panel at identical parameters to read
## against; a structure that is real should survive the sweep, a donor artifact tracks the bird
## whatever the parameters.
##
## Plotted from the embedding by hand rather than through DimPlot(split.by=) for two reasons. Draw
## order matters: whichever bird is drawn last covers the others, and with birds occupying the same
## cell types that is the difference between "bird3 is everywhere" and the truth, so the points are
## shuffled before drawing. And the split panels need the *other* birds' cells behind them in grey,
## without which each panel is a differently-shaped cloud and nothing is comparable across panels.

individual_levels = levels(objs[[1]]$individual)
individual_colors = setNames(
  DiscretePalette_scCustomize(num_colors = length(individual_levels), palette = "varibow"),
  individual_levels
)

umap_df = function(obj, reduction, vars) {
  emb = Embeddings(obj, reduction)[, 1:2]
  colnames(emb) = c("UMAP_1", "UMAP_2")
  bind_cols(as_tibble(emb, rownames = "cell"),
            obj@meta.data[, vars, drop = FALSE] %>% as_tibble())
}

## theme_cowplot leaves plot.background transparent, which is invisible against a white viewer and
## black against a dark one -- fine for the solid-colour panels above, not for these, where the
## point of the figure is grey cells against the background. Paint it white explicitly.
umap_theme_white = function() {
  umap_theme() + theme(plot.background = element_rect(fill = "white", colour = NA))
}

plot_by_individual = function(df, pt.size = 0.4) {
  ggplot(df[sample(nrow(df)), ], aes(UMAP_1, UMAP_2, colour = individual)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_colour_manual(values = individual_colors, drop = FALSE, name = NULL) +
    coord_equal() +
    umap_theme_white() +
    guides(colour = guide_legend(override.aes = list(size = 3)))
}

split_by_individual = function(df, ncol = 3, pt.size = 0.4) {
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

iwalk(objs, function(obj_int_filt, region_cur) {
  reductions = grep("dims", Reductions(obj_int_filt), value = T)

  for (reduction.name in reductions) {
    ind_df = umap_df(obj_int_filt, reduction.name,
                     c("individual", "soup_batch", "soup_library", "position", res_to_use)) %>%
      mutate(individual = droplevels(individual))

    gg = plot_by_individual(ind_df)
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_individual.pdf", region_cur, reduction.name)), gg, base_height = 7, base_asp = 1.15)
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_individual.png", region_cur, reduction.name)), gg, base_height = 7, base_asp = 1.15)

    gg = split_by_individual(ind_df, ncol = 3)
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_individual_split.pdf", region_cur, reduction.name)), gg, base_height = 7, base_asp = 1.5)
    save_plot(file.path(out_dir, sprintf("umap_%s_%s_individual_split.png", region_cur, reduction.name)), gg, base_height = 7, base_asp = 1.5)

    ## Per batch, where there is more than one. Birds are only comparable within a batch and the two
    ## batches cover different regions, so the combined split panels put a two-region cloud next to
    ## another; these restrict the background to the batch's own cells, which is the comparison that
    ## is actually meaningful. png only -- screening panels, not figure panels.
    batches = levels(droplevels(ind_df$soup_batch))
    if (length(batches) > 1) {
      for (b in batches) {
        df_b = ind_df %>% filter(soup_batch == b) %>% mutate(individual = droplevels(individual))
        gg = split_by_individual(df_b, ncol = 3)
        save_plot(file.path(out_dir, sprintf("umap_%s_%s_individual_split_%s.png", region_cur, reduction.name, b)), gg, base_height = 4, base_asp = 2.6)
      }
    }
  }
})

## Cells per bird per cell type -- the count behind every panel above, and the first thing to check
## when a cluster looks donor-driven: a cluster that is one bird's cells says so here before any
## embedding is read.
imap(objs, function(obj_int_filt, region_cur) {
  obj_int_filt@meta.data %>%
    as_tibble() %>%
    count(soup_batch, soup_library, individual, .data[[res_to_use]], .drop = TRUE) %>%
    mutate(region = region_cur, .before = 1)
}) %>%
  bind_rows() %>%
  filter(n > 0) %>%
  write_csv(file.path(out_dir, "cells_per_individual_celltype.csv"))

## The same counts as a per-cluster share of the largest contributing bird -- one number per
## cluster, within batch, for the question the panels are being read for.
cell_shares = map(objs, function(obj_int_filt) {
  obj_int_filt@meta.data %>%
    as_tibble() %>%
    count(soup_batch, .data[[res_to_use]], individual, .drop = TRUE) %>%
    filter(n > 0) %>%
    group_by(soup_batch, .data[[res_to_use]]) %>%
    summarise(n_cells = sum(n),
              n_birds = n_distinct(individual),
              top_bird = individual[which.max(n)],
              top_bird_frac = max(n) / sum(n),
              .groups = "drop")
}) %>%
  bind_rows(.id = "region") %>%
  arrange(-top_bird_frac)
print(as.data.frame(cell_shares), digits = 3)
write_csv(cell_shares, file.path(out_dir, "individual_share_by_celltype.csv"))

# Is a cluster split by bird? ---------------------------------------------

## The share table above answers "is this cluster one bird's cells", which is not the same question
## as "is this cluster split by bird" -- a cluster can draw evenly from three birds and still sit in
## the embedding as one island per bird. Glut-CACNA1H-RA is exactly that case, so the composition
## counts alone would be read the wrong way round.
##
## This is the within-cluster version: of the spread of a cluster's cells in the embedding, how much
## lies between birds rather than within them. It is a plain one-way variance ratio on the two UMAP
## coordinates (which share units, so they are summed unscaled) -- 0 means the birds are
## interleaved, 1 means each bird sits in its own spot. Computed within batch, because birds from
## different batches share no library and the comparison across them is not meaningful, and over
## every reduction, because a split that is real should not depend on the UMAP parameters.
##
## It is descriptive, not a test. Cells are not independent samples: there are three birds per
## batch, one library per dissection, so this ranks clusters for inspection and nothing more.
bird_variance_ratio = function(coords, bird) {
  ss_tot = sum(sweep(coords, 2, colMeans(coords))^2)
  ss_within = sum(vapply(split(seq_len(nrow(coords)), bird), function(i) {
    if (length(i) < 2) return(0)
    sum(sweep(coords[i, , drop = FALSE], 2, colMeans(coords[i, , drop = FALSE]))^2)
  }, numeric(1)))
  if (ss_tot == 0) return(NA_real_)
  1 - ss_within / ss_tot
}

min_cells = 20
bird_split = imap(objs, function(obj_int_filt, region_cur) {
  reductions = grep("dims", Reductions(obj_int_filt), value = T)
  md = obj_int_filt@meta.data

  map(reductions, function(reduction.name) {
    emb = Embeddings(obj_int_filt, reduction.name)[, 1:2]

    md %>%
      as_tibble(rownames = "cell") %>%
      mutate(row = row_number()) %>%
      group_by(soup_batch, .data[[res_to_use]]) %>%
      filter(n() >= min_cells, n_distinct(individual) >= 2) %>%
      summarise(n_cells = n(),
                n_birds = n_distinct(individual),
                bird_var_ratio = bird_variance_ratio(emb[row, , drop = FALSE],
                                                     droplevels(individual)),
                .groups = "drop") %>%
      mutate(reduction = reduction.name, .before = 1)
  }) %>% bind_rows() %>% mutate(region = region_cur, .before = 1)
}) %>% bind_rows()

write_csv(bird_split, file.path(out_dir, "individual_split_by_celltype.csv"))

## Averaged over the sweep, so the ranking is not one parameter choice's opinion.
bird_split_summary = bird_split %>%
  group_by(region, soup_batch, .data[[res_to_use]]) %>%
  summarise(n_cells = max(n_cells),
            n_birds = max(n_birds),
            bird_var_ratio_mean = mean(bird_var_ratio),
            bird_var_ratio_min = min(bird_var_ratio),
            bird_var_ratio_max = max(bird_var_ratio),
            .groups = "drop") %>%
  arrange(-bird_var_ratio_mean)
print(as.data.frame(bird_split_summary), digits = 3)
write_csv(bird_split_summary, file.path(out_dir, "individual_split_by_celltype_summary.csv"))
