#!/usr/bin/env Rscript
# Song-pair DARs with ArchR's bias-matched test, on the consensus regions of the current cisTopic object.
#
# Why not pycisTopic's find_diff_features: it ranks regions by the ratio of group means of *imputed* accessibility, which is a
# linear function of the topic weights, so any topic that differs a little between two groups makes every region loaded on it a
# "DAR" however small the topic is (see the 2026-10 diagnosis: only 20 % of the HVCra-vs-HVCx DARs were supported by the raw
# fragments). ArchR's getMarkerFeatures tests the observed peak counts, Wilcoxon on a background matched to the foreground
# cells for TSSEnrichment and log10(nFrags), the same call multiome/archr/differential_accessibility.qmd uses.
#
# The ArchR project holds its own peak set, so this works on a COPY of it (the original is never written to): the consensus
# regions are added as the copy's peak set and the peak matrix is rebuilt. Steps are skipped when their output exists.
#
#   Rscript make_song_pair_dars_archr.R            # all steps
#
# Outputs (not tracked): the project copy and region_sets_k40_archr/ next to region_sets_k40/ (its DARs_song-pairs/ is the new
# one). song_pair_dars_archr_summary.csv next to this script is tracked.
suppressMessages({library(ArchR); library(GenomicRanges); library(here)})
setwd(here::here())
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))

STORE <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid"
CONSENSUS <- file.path(STORE, "pycisTopic/scATAC/consensus_peak_calling/consensus_regions.bed")
PROJ_COPY <- file.path(STORE, "archr_consensus")
SETS_IN <- file.path(STORE, "pycisTopic/region_sets_k40")
SETS_OUT <- file.path(STORE, "pycisTopic/region_sets_k40_archr")
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
## foreground is the first element; only regions MORE accessible in the foreground are kept, as in the pycisTopic version
## each song cell type against its non-song counterpart: RA vs CACNA1H-1, HVCra vs DACH2-1, HVCx vs DACH2-4
CONTRASTS <- list(c("Glut-CACNA1H-RA", "Glut-CACNA1H-1"), c("Glut-DACH2-HVCra", "Glut-DACH2-1"), c("Glut-DACH2-HVCx", "Glut-DACH2-4"))
CUTOFF <- "FDR <= 0.05 & Log2FC >= 0.585"   # the pycisTopic run used adjusted p <= 0.05 and log2FC >= log2(1.5)

addArchRThreads(threads = THREADS)

## 1. copy of the project (arrow files included); the original is only read
if (!file.exists(file.path(PROJ_COPY, "Save-ArchR-Project.rds"))) {
  src <- loadArchRProject(path.expand(ARCHR_PROJ_DIR), showLogo = FALSE)
  message("copying ", nCells(src), " cells to ", PROJ_COPY)
  saveArchRProject(src, outputDirectory = PROJ_COPY, overwrite = FALSE, load = FALSE)
  rm(src)
}
proj <- loadArchRProject(PROJ_COPY, showLogo = FALSE)
stopifnot(normalizePath(PROJ_COPY) != normalizePath(path.expand(ARCHR_PROJ_DIR)))

## 2. consensus regions as the peak set, then the peak matrix
if (!"PeakMatrix" %in% getAvailableMatrices(proj) || !file.exists(file.path(PROJ_COPY, "consensus_peakset_done"))) {
  b <- read.delim(CONSENSUS, header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(b$V1, IRanges(b$V2 + 1L, b$V3))          # BED is 0-based; region names (chr:start-end) use the BED coordinates
  keep <- as.character(seqnames(gr)) %in% names(seqlengths(getGenomeAnnotation(proj)$chromSizes))
  message("consensus regions: ", length(gr), "; on ArchR chromosomes: ", sum(keep))
  gr <- sort(gr[keep])
  mcols(gr)$name <- paste0(seqnames(gr), ":", start(gr) - 1L, "-", end(gr))
  proj <- addPeakSet(proj, peakSet = gr, force = TRUE)
  proj <- addPeakMatrix(proj, force = TRUE)
  saveArchRProject(proj, load = FALSE)   # persist the new peak set; the markers read the matrix's own coordinates either way
  file.create(file.path(PROJ_COPY, "consensus_peakset_done"))
}
proj <- loadArchRProject(PROJ_COPY, showLogo = FALSE)
message("peak set: ", length(getPeakSet(proj)), " regions")

## 3. hybrid labels (in memory) and the three contrasts
proj <- add_cluster_hybrid(proj)
proj$cluster_hybrid <- as.character(proj$cluster_hybrid)   # getMarkerFeatures errors on a factor
dir.create(file.path(SETS_OUT, "DARs_song-pairs"), recursive = TRUE, showWarnings = FALSE)
## only the current contrasts may be in the folder: SCENIC+ runs motif enrichment on every .bed in it
unlink(list.files(file.path(SETS_OUT, "DARs_song-pairs"), pattern = "\\.bed$", full.names = TRUE))
summ <- list()
for (p in CONTRASTS) {
  nm <- paste0(p[1], "_VS_", p[2])
  rds <- file.path(PROJ_COPY, paste0("markers_", nm, ".rds"))
  if (file.exists(rds)) se <- readRDS(rds) else {
    se <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "cluster_hybrid", useGroups = p[1], bgdGroups = p[2],
                            bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
    saveRDS(se, rds)
  }
  gr <- getMarkers(se, cutOff = CUTOFF, returnGR = TRUE)[[p[1]]]
  bed <- data.frame(chr = as.character(seqnames(gr)), start = start(gr) - 1L, end = end(gr))
  bed <- bed[order(bed$chr, bed$start), ]
  write.table(bed, file.path(SETS_OUT, "DARs_song-pairs", paste0(nm, ".bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  summ[[nm]] <- data.frame(contrast = nm, n_fg_cells = sum(proj$cluster_hybrid == p[1]), n_bg_cells = sum(proj$cluster_hybrid == p[2]),
                           cutoff = CUTOFF, n_dars = nrow(bed))
  message(nm, ": ", nrow(bed), " DARs")
}
write.csv(do.call(rbind, summ), here::here("grn/ra-arco-hvc-nc/pycisTopic/song_pair_dars_archr_summary.csv"), row.names = FALSE)

## 4. region_sets_k40_archr = region_sets_k40 (Topics, DARs_all) + the new DARs_song-pairs
for (d in c("Topics_otsu", "Topics_top_3k", "DARs_all")) {
  if (!dir.exists(file.path(SETS_OUT, d))) file.copy(file.path(SETS_IN, d), SETS_OUT, recursive = TRUE)
}
message("done: ", SETS_OUT)
