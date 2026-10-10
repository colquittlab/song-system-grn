#!/usr/bin/env Rscript
# ArchR DARs for the interneuron side of the MAFB comparison: PVALB-1 and PVALB-2 against the MAFB-low MGE type LAMP5, on the consensus regions
# of the current cisTopic object. Same test and cutoffs as make_song_pair_dars_archr.R (RA vs C1H-1 is made there): bias-matched Wilcoxon on the
# peak matrix of the project copy (archr_consensus), FDR <= 0.05 and Log2FC >= 0.585, regions MORE accessible in the foreground only.
#
# Output goes to dars_extra/, NOT into a region-set folder: SCENIC+ runs motif enrichment on every .bed in a region-set folder.
#
#   Rscript make_interneuron_dars_archr.R
suppressMessages({library(ArchR); library(GenomicRanges); library(here)})
setwd(here::here())
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))
STORE <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid"
PROJ_COPY <- file.path(STORE, "archr_consensus")
OUT <- file.path(STORE, "pycisTopic/dars_extra")
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
CUTOFF <- "FDR <= 0.05 & Log2FC >= 0.585"
CONTRASTS <- list(c("GABA-MGE-PVALB-1", "GABA-MGE-LAMP5"), c("GABA-MGE-PVALB-2", "GABA-MGE-LAMP5"))
addArchRThreads(threads = THREADS)
stopifnot(file.exists(file.path(PROJ_COPY, "Save-ArchR-Project.rds")))
proj <- loadArchRProject(PROJ_COPY, showLogo = FALSE)
proj <- add_cluster_hybrid(proj)
proj$cluster_hybrid <- as.character(proj$cluster_hybrid)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
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
  write.table(bed, file.path(OUT, paste0(nm, ".bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  summ[[nm]] <- data.frame(contrast = nm, n_fg_cells = sum(proj$cluster_hybrid == p[1]), n_bg_cells = sum(proj$cluster_hybrid == p[2]), cutoff = CUTOFF, n_dars = nrow(bed))
  message(nm, ": ", nrow(bed), " DARs")
}
write.csv(do.call(rbind, summ), here::here("grn/ra-arco-hvc-nc/pycisTopic/interneuron_dars_archr_summary.csv"), row.names = FALSE)
