## Differential accessibility GABA-MGE-PVALB-1 vs GABA-LGE-1 (old names GABA-3 vs GABA-1-1), for the MAFB-domain
## element analysis (comparative_genomics/phyloP_matched/loci_focus), which used PVALB-2 (GABA-4-1) until Fig. 6 was
## aligned with Fig. 4 (PVALB-1 is the GABAergic comparator there).
##
## The test is the one peaks_gaba.qmd runs for each pair in hybrid_pairs_gaba (PeakMatrix, cluster_hybrid, Wilcoxon,
## bias TSSEnrichment + log10(nFrags)), on the project's fixed peak set; the table has the columns of the existing
## panelE_inputs/dar_*.tsv (export_dar.R): every peak, no cutoffs. The SnapATAC2 MACS3 peaks for GABA-3 are already
## called for every cluster in snapatac/call_peaks.ipynb; this copies that file next to the table.
## Writes, to multiome/archr/peaks_gaba_hybrid/:
##   dar_pvalb1_lge1.tsv        chrom start end Log2FC FDR, every peak
##   dar_pvalb1_lge1_se.rds     the getMarkerFeatures SummarizedExperiment
##   snapatac_GABA-3.bed        SnapATAC2 MACS3 peaks for GABA-3 (= GABA-MGE-PVALB-1)
##
##   Rscript multiome/archr/dar_pvalb1_lge1.R        (on lark; ~15-30 min, 36 threads)
suppressMessages({library(ArchR); library(tidyverse); library(SummarizedExperiment)})
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))

SNAP_PEAKS <- "/hdd/jupyter/brad/snapatac/motor-pathway_multiome/motor-pathway_multiome_cellbender_cluster-snrna-cr/call_peaks"
out_dir <- here::here("multiome/archr", "peaks_gaba_hybrid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
addArchRThreads(threads = 36)

proj <- add_cluster_hybrid(loadArchRProject(ARCHR_PROJ_DIR))
pair <- list("GABA-MGE-PVALB-1", "GABA-LGE-1")
stopifnot(all(unlist(pair) %in% unique(proj$cluster_hybrid)), "PeakMatrix" %in% getAvailableMatrices(proj))
message(sprintf("[%s] cells: %s", Sys.time(), paste(sprintf("%s %d", unlist(pair),
        sapply(unlist(pair), function(g) sum(proj$cluster_hybrid == g))), collapse = ", ")))

se <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "cluster_hybrid",
                        useGroups = pair[[1]], bgdGroups = pair[[2]],
                        bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
saveRDS(se, file.path(out_dir, "dar_pvalb1_lge1_se.rds"))

rd <- rowData(se)
df <- data.frame(chrom = rd$seqnames, start = rd$start, end = rd$end,
                 Log2FC = assay(se, "Log2FC")[, 1], FDR = assay(se, "FDR")[, 1])
write.table(df, file.path(out_dir, "dar_pvalb1_lge1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dom <- df$chrom == "chr20" & df$end > 12551036 & df$start < 13549500
message(sprintf("  %d peaks; FDR < 0.1 & Log2FC > 0: %d genome-wide, %d in the MAFB domain (chr20:12,551,036-13,549,500)",
                nrow(df), sum(df$FDR < 0.1 & df$Log2FC > 0, na.rm = TRUE),
                sum(dom & df$FDR < 0.1 & df$Log2FC > 0, na.rm = TRUE)))

src <- file.path(SNAP_PEAKS, "GABA-3.bed")
if (file.exists(src)) {
  file.copy(src, file.path(out_dir, "snapatac_GABA-3.bed"), overwrite = TRUE)
  message(sprintf("  copied %s (%d peaks)", src, length(readLines(src))))
} else {
  message(sprintf("  NOT FOUND: %s -- available: %s", src, paste(list.files(SNAP_PEAKS, "\\.bed$"), collapse = ", ")))
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "dar_pvalb1_lge1_sessionInfo.txt"))
