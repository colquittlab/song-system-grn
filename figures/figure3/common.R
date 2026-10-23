# figures/figure3/common.R -- inputs and shared helpers for the Figure 3 panel scripts.
#
# Every panel script sources this file (and the repository's config/ files through it) and nothing else in this folder, so each panel can be
# regenerated on its own. Inputs are the outputs of the upstream pipelines named in README.md; edit the paths below (or set the environment variables)
# to point at them. Panels are written to FIG3_OUT (default: figures/figure3/output/, not tracked).
suppressMessages({library(here); library(tidyverse)})
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))

FIG3_OUT <- Sys.getenv("FIG3_OUT", here::here("figures/figure3/output"))
dir.create(FIG3_OUT, showWarnings = FALSE, recursive = TRUE)

# ---- inputs -----------------------------------------------------------------------------------------------------------------------------------------
# Hybrid-labeled multiome Seurat object with the UMAP reductions (multiome/seurat/reduction_viz/combined_all_umap_hybrid.R)
SEURAT_HYBRID_OBJ <- Sys.getenv("FIG3_SEURAT_OBJ", here::here("multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"))

# ---- cell types -------------------------------------------------------------------------------------------------------------------------------------
# The three song-system projection neurons and their matched non-song neurons, in the row order of the panel C genome-browser tracks; colors are the
# track colors of the browser export
SONG_PAIRS <- tribble(~cluster, ~label, ~color,
                      "Glut-CACNA1H-RA",  "RA",    "#e31a1c",
                      "Glut-CACNA1H-1",   "C1H-1", "#fb9a99",
                      "Glut-DACH2-HVCra", "HVCra", "#ff7f00",
                      "Glut-DACH2-1",     "D2-1",  "#fdbf6f",
                      "Glut-DACH2-HVCx",  "HVCx",  "#ff7f00",
                      "Glut-DACH2-4",     "D2-4",  "#fdbf6f")

# ---- expression ------------------------------------------------------------------------------------------------------------------------------------
#' SCT log-normalized expression of one gene per cell type: mean of per-bird means (a bird = one souporcell assignment within one library, kept when it
#' has >= min_cells cells of that type), with the spread across birds and the per-cell fraction expressing. Returns list(summary, per_bird).
gene_expression_by_bird <- function(obj, gene, clusters, min_cells = 10) {
  suppressMessages(library(Seurat))
  DefaultAssay(obj) <- "SCT"; stopifnot(gene %in% rownames(obj))
  d <- FetchData(obj, vars = c(gene, "cluster_hybrid", "sample_short", "assignment"), layer = "data") %>% rename(expr = 1) %>% filter(cluster_hybrid %in% clusters)
  bird <- d %>% group_by(cluster = cluster_hybrid, sample_short, assignment) %>% summarize(n_cells = n(), pseudobulk = mean(expr), .groups = "drop") %>%
    filter(n_cells >= min_cells)
  summ <- tibble(cluster = clusters) %>%
    left_join(bird %>% group_by(cluster) %>% summarize(n_birds = n(), mean_expr = mean(pseudobulk), sd_birds = sd(pseudobulk), .groups = "drop"), by = "cluster") %>%
    left_join(d %>% group_by(cluster = cluster_hybrid) %>% summarize(n_cells = n(), mean_expr_cells = mean(expr), pct_expressing = 100 * mean(expr > 0), .groups = "drop"), by = "cluster")
  list(summary = summ, per_bird = bird)
}
