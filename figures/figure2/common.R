# figures/figure2/common.R -- inputs and shared helpers for the Figure 2 R panel scripts (A, B, C).
#
# Every R panel script sources this file (and the repository's config/ files through it) and nothing else in this folder, so each panel can be
# regenerated on its own. Inputs are the outputs of the upstream pipelines named in README.md; edit the paths below (or set the environment variables)
# to point at them. Panels are written to FIG2_OUT (default: figures/figure2/output/, not tracked).
#
# The dissection-position colours the original scripts took from the lab's shared utilities (/opt/colquittlab/utils/R/common_aesthetics.R) are
# copied below, so these scripts do not need that directory.
suppressMessages({library(here); library(tidyverse)})
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))
select <- dplyr::select

FIG2_OUT <- Sys.getenv("FIG2_OUT", here::here("figures/figure2/output"))
dir.create(FIG2_OUT, showWarnings = FALSE, recursive = TRUE)

# ---- inputs -----------------------------------------------------------------------------------------------------------------------------------------
# snRNA-seq object of the glutamatergic cells with the hybrid labels (snrna/trees/celltypes_hclust_glut_robustness_hybrid.qmd)
GLUT_OBJ <- Sys.getenv("FIG2_GLUT_OBJ", here::here("snrna/trees/celltypes_hclust_glut_robustness_hybrid/obj_clustered_glutamatergic.qs2"))
# RCTD label transfer from the snRNA-seq hybrid labels onto the proseg cells (xenium/label_transfer/hpc_rctd_proseg_hybrid/): per-cell calls
# (rctd_all.csv.gz), centroids (proseg_cell_metadata.csv.gz) and transcript counts (proseg_ncount.csv.gz)
RCTD_DIR <- Sys.getenv("FIG2_RCTD_DIR", file.path(path.expand(XENIUM_PROC_DIR), "hpc_rctd_proseg_hybrid"))

# ---- dissection positions ---------------------------------------------------------------------------------------------------------------------------
POSITION_COLORS <- c(hvc = "#ff7f00", nc = "#fdbf6f", lman = "#1f78b4", nr = "#a6cee3", ra = "#e31a1c", arco = "#fb9a99")

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

#' Every proseg cell with its section, centroid, RCTD call (spot_class, first_type) and transcript count (nCount)
read_xenium_calls <- function(dir = RCTD_DIR) {
  md <- read_csv(file.path(dir, "proseg_cell_metadata.csv.gz"), show_col_types = FALSE) %>% select(cell, section_id, x_centroid, y_centroid)
  rctd <- read_csv(file.path(dir, "rctd_all.csv.gz"), show_col_types = FALSE) %>% select(cell, spot_class, first_type)
  nc <- read_csv(file.path(dir, "proseg_ncount.csv.gz"), show_col_types = FALSE)
  md %>% left_join(rctd, by = "cell") %>% left_join(nc, by = "cell")
}
