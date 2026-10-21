# figures/figure1/common.R -- inputs and shared helpers for the Figure 1 panel scripts.
#
# Every panel script sources this file (and the repository's config/ files through it) and nothing else in this folder, so each panel can be
# regenerated on its own. Inputs are the outputs of the upstream pipelines named in README.md; edit the paths below (or set the environment variables)
# to point at them. Panels are written to FIG1_OUT (default: figures/figure1/output/, not tracked).
#
# The two pieces the original analysis scripts took from the lab's shared utilities (/opt/colquittlab/utils/R/) -- calc_tissue_specificity() from
# stats.R and the dissection-position colours from common_aesthetics.R -- are copied below, so these scripts do not need that directory.
suppressMessages({library(here); library(tidyverse)})
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))
select <- dplyr::select

FIG1_OUT <- Sys.getenv("FIG1_OUT", here::here("figures/figure1/output"))
dir.create(FIG1_OUT, showWarnings = FALSE, recursive = TRUE)

# ---- inputs -----------------------------------------------------------------------------------------------------------------------------------------
# snRNA-seq Seurat object with the hybrid division-based cell-type labels (celltype_hybrid; snrna/naming/hybrid_division_naming.qmd)
SNRNA_HYBRID_OBJ <- Sys.getenv("FIG1_SNRNA_OBJ", here::here("snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2"))
# Proseg segmentation of the ten Xenium sections, one folder per section (xenium/preprocessing/hpc_proseg/)
PROSEG_DIR <- Sys.getenv("FIG1_PROSEG_DIR", path.expand(XENIUM_PROSEG_DIR))
# Transcripts binned at 10 um per section (panel_C_bin_transcripts.py writes tx_density_<section>_10um.csv.gz here if it is missing)
TX_DENSITY_DIR <- Sys.getenv("FIG1_TX_DENSITY_DIR", here::here("xenium/transcript_density"))
# RCTD label transfer from the snRNA-seq hybrid labels onto the proseg cells (xenium/label_transfer/hpc_rctd_proseg_hybrid/): per-cell calls
# (rctd_all.csv.gz), centroids (proseg_cell_metadata.csv.gz) and transcript counts (proseg_ncount.csv.gz)
RCTD_DIR <- Sys.getenv("FIG1_RCTD_DIR", file.path(path.expand(XENIUM_PROC_DIR), "hpc_rctd_proseg_hybrid"))

# ---- cell types -------------------------------------------------------------------------------------------------------------------------------------
CELLTYPE_COL <- "celltype_hybrid"
# Row order of panel G (snrna/trees/celltypes_hclust_all_hybrid.qmd, snrna/deg/dotplot_manual_markers_hybrid.qmd); the 47 labelled types
CT_ORDER <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCra-Int", "Glut-DACH2-HVCx",
              "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
              "Glut-DACH2-LMANsh", "Glut-DACH2-LMANco",
              "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8", "Glut-SATB2-1",
              "Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
              "Glut-GABA",
              "Glut-Im", "Glut-NB", "Glut-NSC",
              "GABA-LGE-1", "GABA-LGE-2",
              "GABA-MGE-SST-1", "GABA-MGE-SST-2", "GABA-MGE-SST-3",
              "GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2", "GABA-MGE-PVALB-3",
              "GABA-MGE-LAMP5", "GABA-MGE-ST18", "GABA-MGE-LHX8",
              "GABA-CGE-1", "GABA-CGE-2",
              "GABA-Im",
              "Astro-1", "Astro-2", "Astro-3", "Epen",
              "Oligo-1", "Oligo-2", "Oligo-3", "OPC",
              "Micro", "Endo")
# Dissection positions, in panel G's column order
POSITIONS <- c("hvc", "nc", "lman", "nr", "ra", "arco")
POSITION_COLORS <- c(hvc = "#ff7f00", nc = "#fdbf6f", lman = "#1f78b4", nr = "#a6cee3", ra = "#e31a1c", arco = "#fb9a99")

# Bar colours of panel G: glutamatergic / glutamatergic precursor / GABAergic / GABAergic precursor / non-neuronal. Glut-GABA counts as glutamatergic.
CLASS_COLORS <- c("Glut" = "#008000", "Glut-Pre" = "#90ca90", "GABA" = "#800000", "GABA-Pre" = "#c37474", "Non-neuron" = "#808000")
celltype_class <- function(ct) case_when(ct %in% c("Glut-NSC", "Glut-NB", "Glut-Im") ~ "Glut-Pre",
                                         grepl("^Glut", ct) ~ "Glut",
                                         ct %in% "GABA-Im" ~ "GABA-Pre",
                                         grepl("^GABA", ct) ~ "GABA",
                                         TRUE ~ "Non-neuron")

#' snRNA-seq object restricted to the labelled cells (celltype_hybrid is NA for the two clusters excluded as artifacts, Glut-Arco-4 and Glut-Nido-3)
read_snrna <- function(fname = SNRNA_HYBRID_OBJ) {
  suppressMessages({library(Seurat); library(qs2)})
  obj <- qs_read(fname, nthreads = 8)
  obj <- subset(obj, cells = colnames(obj)[!is.na(obj@meta.data[[CELLTYPE_COL]])])
  stopifnot("cell types differ from CT_ORDER" = setequal(CT_ORDER, unique(obj@meta.data[[CELLTYPE_COL]])))
  DefaultAssay(obj) <- "SCT"
  obj
}

#' Tissue specificity of a vector of proportions: KL divergence (nats) from the uniform distribution (/opt/colquittlab/utils/R/stats.R)
calc_tissue_specificity <- function(x, scale = FALSE) {
  if (scale) x <- (x - min(x)) / (max(x) - min(x))
  x_norm <- x / sum(x)
  y <- x_norm * log(x_norm / (1 / length(x)))
  y[is.nan(y)] <- 0
  sum(y)
}

# ---- Xenium -----------------------------------------------------------------------------------------------------------------------------------------
SECTIONS <- c("OR52YW26_1_4", "OR52YW26_1_7", "OR52YW26_2_2", "OR52YW26_2_4", "OR52YW26_2_7",
              "OR69PU4_1_4", "OR69PU4_1_7", "OR69PU4_2_2", "OR69PU4_2_4", "OR69PU4_2_7")
# Per-section rotation/flip to a common orientation (posterior left, dorsal up). Rotations in the standard sense (x right, y up): ccw90 is
# (x, y) -> (-y, x), cw90 is (x, y) -> (y, -x); a flip negates x after any rotation.
SECTION_TRANSFORM <- tribble(
  ~section_id,     ~rotate,  ~flip_h,
  "OR52YW26_1_4",  "none",   FALSE,
  "OR52YW26_1_7",  "ccw90",  FALSE,
  "OR52YW26_2_2",  "ccw90",  TRUE,
  "OR52YW26_2_4",  "ccw90",  FALSE,
  "OR52YW26_2_7",  "cw90",   TRUE,
  "OR69PU4_1_4",   "none",   FALSE,
  "OR69PU4_1_7",   "ccw90",  TRUE,
  "OR69PU4_2_2",   "cw90",   FALSE,
  "OR69PU4_2_4",   "ccw90",  FALSE,
  "OR69PU4_2_7",   "ccw90",  FALSE)

#' Orient x/y of one section's table
orient_xy <- function(df, section, xcol = "x", ycol = "y") {
  tf <- SECTION_TRANSFORM[SECTION_TRANSFORM$section_id == section, ]
  stopifnot("section not in SECTION_TRANSFORM" = nrow(tf) == 1)
  x <- df[[xcol]]; y <- df[[ycol]]
  xr <- switch(tf$rotate, none = x, ccw90 = -y, cw90 = y)
  yr <- switch(tf$rotate, none = y, ccw90 = x, cw90 = -x)
  df[[xcol]] <- if (tf$flip_h) -xr else xr
  df[[ycol]] <- yr
  df
}

# RCTD calls counted as confident: singlets and confident doublets
CONFIDENT_CLASSES <- c("singlet", "doublet_certain")
# Transcript floor for the cells shown in panels F and G (composition below ~400 transcripts is increasingly depth-driven). The figure uses 400;
# set FIG1_MIN_TRANSCRIPTS=0 for the version with every confident call. Outputs that depend on it are tagged _gte<N> or _allcells.
MIN_TRANSCRIPTS <- as.numeric(Sys.getenv("FIG1_MIN_TRANSCRIPTS", "400"))
DEPTH_TAG <- if (MIN_TRANSCRIPTS > 0) paste0("_gte", MIN_TRANSCRIPTS) else "_allcells"

#' Every proseg cell with its section, centroid, RCTD call (spot_class, first_type) and transcript count (nCount)
read_xenium_calls <- function(dir = RCTD_DIR) {
  md <- read_csv(file.path(dir, "proseg_cell_metadata.csv.gz"), show_col_types = FALSE) %>% select(cell, section_id, x_centroid, y_centroid)
  rctd <- read_csv(file.path(dir, "rctd_all.csv.gz"), show_col_types = FALSE) %>% select(cell, spot_class, first_type)
  nc <- read_csv(file.path(dir, "proseg_ncount.csv.gz"), show_col_types = FALSE)
  md %>% left_join(rctd, by = "cell") %>% left_join(nc, by = "cell")
}
