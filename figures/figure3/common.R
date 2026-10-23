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

# ---- ArchR (panels D, E) ------------------------------------------------------------------------------------------------------------------------------
# Panel scripts attach ArchR themselves; it attaches plyr after the tidyverse, so they rebind the masked dplyr verbs after loading it.
# list(ra, hvc_1, hvc_2) of getMarkerFeatures() results, saved by multiome/archr/peaks_glut.qmd ("Find DARs"). When the file is absent the three tests
# are recomputed exactly as there and cached in FIG3_CACHE (song_pair_markers()).
DAR_MARKERS <- Sys.getenv("FIG3_DAR_MARKERS", here::here("multiome/archr/peaks_glut_hybrid/peaks_nonsong_glut_se.rds"))
# scratch for the values computed from the shared ArchR project (the project itself is only read, never saved)
FIG3_CACHE <- Sys.getenv("FIG3_CACHE", file.path(FIG3_OUT, "cache"))
dir.create(FIG3_CACHE, showWarnings = FALSE, recursive = TRUE)

# ArchR writes its logs (ArchRLogs/) and temporary files under the working directory: run its calls from the cache
in_cache <- function(expr) { wd <- setwd(FIG3_CACHE); on.exit(setwd(wd)); expr }
# The shared project with the hybrid labels added in memory (multiome/archr/hybrid_labels.R), loaded only when a cached value is missing.
# addArchRThreads(threads > 1) also switches R to the L'Ecuyer-CMRG generator, as in the notebook, whose results depend on it.
.proj <- NULL
archr_hybrid <- function() {
  if (is.null(.proj)) {
    source(here::here("multiome/archr/hybrid_labels.R"))
    addArchRThreads(threads = 16)
    .proj <<- in_cache(add_cluster_hybrid(loadArchRProject(ARCHR_PROJ_DIR, showLogo = FALSE)))
  }
  .proj
}
#' The three song vs. non-song marker tests on the PeakMatrix (bias-matched Wilcoxon; ArchR seeds its bias matching with seed 1)
song_pair_markers <- function() {
  if (file.exists(DAR_MARKERS)) return(readRDS(DAR_MARKERS))
  f <- file.path(FIG3_CACHE, "peaks_nonsong_glut_se.rds")
  if (!file.exists(f)) {
    message(DAR_MARKERS, " not found: recomputing the marker tests into ", f)
    proj <- archr_hybrid()
    m <- in_cache(map(hybrid_pairs_glut, \(p) getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "cluster_hybrid", useGroups = p[[1]],
                                                                 bgdGroups = p[[2]], bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")))
    saveRDS(m, f)
  }
  readRDS(f)
}

# ---- song vs. non-song expression of TF families (panels G, H) --------------------------------------------------------------------------------------
# Pseudobulk DESeq2 table of multiome/seurat/deg/song-surround_deg_glut.qmd, and the hybrid-labeled object chromvar.qmd took its gene names from
DEG_TABLE <- Sys.getenv("FIG3_DEG_TABLE", file.path(MULTIOME_SEURAT_DIR, "deg/song-surround_deg_glut/marker_genes_cluster_SCT.rds"))
HYBRID_LABEL_OBJ <- Sys.getenv("FIG3_HYBRID_LABEL_OBJ", here::here("multiome/seurat/label_transfer_hybrid/obj_clustered_hybrid.qs2"))
# TF-family patterns of chromvar.qmd (the list in force when the points_deg and heatmap_deg files were written)
TF_TERMS <- list(AP1 = c("^SMARCC", "^JUN", "^FOS", "^JDP"), NFI = c("^NFI"), IRF = c("^IRF[0-9]+$"), MAF = c("^MAF[^0-9]", "^MAF$"), NR3 = c("^NR3C"),
                 CIC = "^CIC", BACH = "^BACH", BATF = "^BATF", BHLH = "^BHLH", OLIG = "^OLIG", NEURO = "^NEURO", ATOH = "^ATOH", TAL = "^TAL[0-9]+$",
                 MEF2 = "^MEF2", PPAR = "^PPAR", SIX = "^SIX", TCF = c("^TCF[0-9]+", "^TCFL[0-9]+"), AR = "^AR$", PGR = "^PGR$", ARID = "^ARID",
                 EMX = "^EMX", RFX = "^RFX[0-9]+")
#' The DEG table (computed under the pre-hybrid labels), its `compare` column carried over to the hybrid names as in the notebook
read_song_deg <- function() {
  source(here::here("multiome/archr/hybrid_labels.R"))
  deg <- readRDS(DEG_TABLE)
  deg$compare <- unname(hybrid_cluster_map()[as.character(deg$compare)])
  stopifnot(!anyNA(deg$compare))
  deg
}
#' Gene names of the object (Features(obj_glut) in the notebook: the default assay's genes)
hybrid_features <- function() { f <- rownames(qs2::qs_read(HYBRID_LABEL_OBJ, nthreads = 8)); invisible(gc()); f }

# ---- expression (panel C) --------------------------------------------------------------------------------------------------------------------------

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
