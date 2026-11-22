## Peak-to-gene links reaching across the whole MAFB (DHX35-TOP1) domain.
##
## peaks2geneslinks.qmd calls addPeak2GeneLinks with ArchR's default maxDist of 250 kb, so MAFB links stop at
## ~13.12-13.61 Mb, while the domain (chr20:12,551,036-13,547,533) reaches ~820 kb upstream of the MAFB TSS
## (13,369,534): most of its DARs were never tested. This reruns the same call (same project, hybrid labels,
## Harmony_ATAC, GeneExpressionMatrix, seed) with maxDist = 1 Mb for
##   ra           Glut-CACNA1H-RA + Glut-CACNA1H-1 cells   (as in peaks2geneslinks_hybrid)
##   gaba         all GABA clusters                      (as in peaks2geneslinks_hybrid)
##   pvalb2_lge1  GABA-MGE-PVALB-2 + GABA-LGE-1 cells       (the G4 DAR contrast)
## and writes, to multiome/archr/peaks2geneslinks_hybrid_maxdist1Mb/:
##   peak2geneLinks_<group>_corCutOff-<c>.bed     the usual link BEDs (c = 0.2..0.45; FDR <= 0.01, VarQ > 0.25),
##                                                 same format as peaks2geneslinks_hybrid, genome-wide
##   mafb_domain_pairs_<group>.tsv                 EVERY tested peak-gene pair with a peak in the domain, with
##                                                 Correlation, FDR and VarQ -- no cutoffs, so each DAR is assessed
## Correlations are recomputed on fresh low-overlap cell aggregates, so pairs inside 250 kb can differ slightly
## from the existing files; FDR is over the larger (1 Mb) test set.
##
##   Rscript multiome/archr/peaks2geneslinks_mafb_domain.R        (on lark; ~1-2 h per group, 36 threads)
suppressMessages({library(ArchR); library(tidyverse)})
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))

MAXDIST <- 1e6
DOMAIN <- GRanges("chr20", IRanges(12551036, 13547533))
out_dir <- here::here("multiome/archr", "peaks2geneslinks_hybrid_maxdist1Mb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
addArchRThreads(threads = 36)

proj <- add_cluster_hybrid(loadArchRProject(ARCHR_PROJ_DIR))
groups <- list(ra = hybrid_pairs_glut$ra, gaba = hybrid_groups$gaba, pvalb2_lge1 = hybrid_pairs_gaba$pvalb2_lge1)
stopifnot(all(unlist(groups) %in% unique(proj$cluster_hybrid)))
cell_names <- getCellNames(proj)

## as getPeak2GeneLinksAsBed in peaks2geneslinks.qmd
getPeak2GeneLinksAsBed <- function(proj, resolution = 100, corCutOff = .45, FDRCutOff = 1e-02,
                                   varCutOffATAC = 0.25, varCutOffRNA = 0.25) {
  p2g <- metadata(proj@peakSet)$Peak2GeneLinks
  p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), , drop = FALSE]
  if (!is.null(varCutOffATAC)) p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC), ]
  if (!is.null(varCutOffRNA)) p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA), ]
  peakSet <- metadata(p2g)$peakSet
  peak_name <- paste(seqnames(peakSet), paste(start(peakSet), end(peakSet), sep = "-"), sep = "_")
  geneSet <- metadata(p2g)$geneSet
  peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
  geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
  summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
  geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
  loops <- ArchR:::.constructGR(seqnames = seqnames(peakSummits[p2g$idxATAC]),
                                start = summitTiles[p2g$idxATAC], end = geneTiles[p2g$idxRNA])
  mcols(loops)$value <- p2g$Correlation
  mcols(loops)$FDR <- p2g$FDR
  mcols(loops)$geneID <- geneSet$name[p2g$idxRNA]
  mcols(loops)$peakID <- peak_name[p2g$idxATAC]
  loops <- loops[order(mcols(loops)$value, decreasing = TRUE)]
  loops <- unique(loops)
  loops <- loops[width(loops) > 0]
  loops <- sort(sortSeqlevels(loops))
  data.frame(chrom = seqnames(loops), chromStart = start(loops), chromEnd = end(loops),
             name = paste(loops$geneID, loops$peakID, sep = "<>"), score = round(1000 * loops$value),
             value = loops$value, exp = ".", color = "red",
             sourceChrom = seqnames(loops), sourceStart = start(loops) - 50, sourceEnd = start(loops) + 50,
             sourceName = loops$geneID, sourceStrand = ".",
             targetChrom = seqnames(loops), targetStart = end(loops) - 50, targetEnd = end(loops) + 50,
             targetName = loops$peakID, targetStrand = ".")
}

for (g in names(groups)) {
  message(sprintf("[%s] %s: addPeak2GeneLinks, maxDist = %g", Sys.time(), g, MAXDIST))
  p <- addPeak2GeneLinks(ArchRProj = proj, cellsToUse = cell_names[proj$cluster_hybrid %in% unlist(groups[[g]])],
                         reducedDims = "Harmony_ATAC", useMatrix = "GeneExpressionMatrix", maxDist = MAXDIST)
  for (cc in seq(0.2, .45, .05))
    write_delim(getPeak2GeneLinksAsBed(p, corCutOff = cc),
                file.path(out_dir, sprintf("peak2geneLinks_%s_corCutOff-%s.bed", g, cc)), delim = "\t", col_names = FALSE)
  ## every tested pair whose peak lies in the domain, no cutoffs
  p2g <- metadata(p@peakSet)$Peak2GeneLinks
  pk <- metadata(p2g)$peakSet; gs <- metadata(p2g)$geneSet
  in_dom <- overlapsAny(pk[p2g$idxATAC], DOMAIN)
  d <- tibble(peak_chrom = as.character(seqnames(pk))[p2g$idxATAC], peak_start = start(pk)[p2g$idxATAC],
              peak_end = end(pk)[p2g$idxATAC], gene = gs$name[p2g$idxRNA], gene_tss = start(resize(gs, 1, "start"))[p2g$idxRNA],
              correlation = p2g$Correlation, fdr = p2g$FDR, varq_atac = p2g$VarQATAC, varq_rna = p2g$VarQRNA)[in_dom, ]
  write_tsv(d, file.path(out_dir, sprintf("mafb_domain_pairs_%s.tsv", g)))
  m <- filter(d, gene == "MAFB")
  message(sprintf("  %s: %d domain pairs tested; MAFB: %d peaks tested (%.0f-%.0f), %d with r >= 0.2 & FDR <= 0.01",
                  g, nrow(d), nrow(m), min(m$peak_start), max(m$peak_end), sum(m$correlation >= 0.2 & m$fdr <= 0.01)))
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
