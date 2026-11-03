#!/usr/bin/env Rscript
# 01_reference_signatures.R -- adult song-neuron gene signatures for scoring the chicken OE data.
#
# For each axis (song cluster vs its non-song sister; AXES in oe_helpers.R):
#   1. pseudobulk DESeq2 on the adult finch snRNA reference, replicate = library x bird x cell type (>= 10 cells),
#      design ~ lib + grp when estimable, else ~ grp (as in figures/figure4/common.R);
#   2. a specificity gate (songbird_dev nido/ session): a gene up in the song cluster must also separate the song cluster
#      from ALL adult glutamatergic cells (AUROC >= 0.70), a gene down must sit low in it (AUROC <= 0.30). A pairwise
#      contrast alone is not a marker set -- a gene can beat the sister and still be ordinary across the pallium;
#   3. mapping to chicken by OMA 1:1 orthologs (the snrna/integration/ rule).
# Signature sets (up-only scoring): <axis>__SONG_STRICT (up + gate), <axis>__SONG_ALL (up), <axis>__SISTER_STRICT
# (down + gate, i.e. genes the song neuron has LOST relative to its sister), <axis>__SISTER_ALL.
#
#   Rscript overexpression/reanalysis/01_reference_signatures.R
# Writes tables/reference_signature_genes.csv (tracked) and output/reference_de_<axis>.csv.gz (all genes).
suppressMessages({library(Seurat); library(qs2)})
source(here::here("overexpression/reanalysis/oe_helpers.R"))
set.seed(2026)

obj <- qs_read(ref_obj_path(), nthreads = 16)
meta <- obj@meta.data
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
rm(obj); invisible(gc())

keep <- !grepl("/", meta$assignment) & !is.na(meta$celltype_hybrid)
meta <- meta[keep, ]; cnt <- as(cnt[, keep], "dgCMatrix")   # the layer is BPCells-backed on disk
cat(sprintf("reference: %d nuclei, %d genes, %d libraries\n", ncol(cnt), nrow(cnt), n_distinct(meta$sample_name)))

## 1. pseudobulk DESeq2 per axis
f <- paste(meta$sample_name, meta$assignment, meta$celltype_hybrid, sep = "|")
pbs <- pseudobulk_sum(cnt, f)
parts <- do.call(rbind, strsplit(names(pbs$n), "|", fixed = TRUE))
info <- tibble(id = names(pbs$n), lib = parts[, 1], bird = parts[, 2], cluster = parts[, 3], n = pbs$n)

de <- pmap(AXES, function(axis, song, sister) {
  s <- info %>% filter(cluster %in% c(song, sister), n >= MIN_CELLS)
  cd <- data.frame(row.names = s$id, lib = factor(s$lib), grp = factor(ifelse(s$cluster == song, "song", "sister"), levels = c("sister", "song")))
  design <- estimable_design(cd, ~ lib + grp, ~ grp)
  r <- deseq_pb(pbs$pb, cd, design, contrast = c("grp", "song", "sister"))
  cat(sprintf("  %-6s %s vs %s: %d vs %d pseudobulk samples (%d libraries), design %s, %d genes, %d up / %d down (padj < %.2f, |lfc| >= %.2f)\n",
              axis, song, sister, sum(cd$grp == "song"), sum(cd$grp == "sister"), n_distinct(s$lib), format(design), nrow(r),
              sum(r$padj < PADJ_MAX & r$lfc >= LFC_MIN, na.rm = TRUE), sum(r$padj < PADJ_MAX & r$lfc <= -LFC_MIN, na.rm = TRUE),
              PADJ_MAX, LFC_MIN))
  r %>% mutate(axis = axis)
}) %>% set_names(AXES$axis)

## 2. specificity gate across adult glutamatergic cells (log-normalized RNA)
glut <- grepl("^Glut-", meta$celltype_hybrid) & !meta$celltype_hybrid %in% ADULT_GLUT_EXCLUDE
cat(sprintf("adult glutamatergic cells for the gate: %d in %d clusters\n", sum(glut), n_distinct(meta$celltype_hybrid[glut])))
cand <- unique(unlist(map(de, ~ .x$gene[!is.na(.x$padj) & .x$padj < PADJ_MAX & abs(.x$lfc) >= LFC_MIN])))
lc <- cnt[cand, glut]
lc <- log1p(t(t(lc) / colSums(cnt[, glut])) * 1e4)
gl <- meta$celltype_hybrid[glut]

sig <- imap_dfr(de, function(r, axis) {
  song <- AXES$song[AXES$axis == axis]
  pos <- gl == song
  r %>% filter(!is.na(padj), padj < PADJ_MAX, abs(lfc) >= LFC_MIN) %>%
    mutate(adult_auroc = map_dbl(gene, ~ fast_auroc(as.numeric(lc[.x, ]), pos)),
           direction = ifelse(lfc > 0, "SONG", "SISTER"),
           specific = ifelse(lfc > 0, adult_auroc >= MIN_AUROC, adult_auroc <= 1 - MIN_AUROC))
})

## 3. chicken orthologs
ck2zf <- oma_ck_to_zf()
sig <- sig %>% mutate(chicken = names(ck2zf)[match(gene, ck2zf)], has_ortholog = !is.na(chicken),
                      excluded_tf = gene %in% TRANSGENE_TFS)

summ <- sig %>% group_by(axis, direction) %>%
  summarise(de = n(), specific_ortholog = sum(specific & has_ortholog & !excluded_tf), ortholog = sum(has_ortholog & !excluded_tf),
            specific = sum(specific), .groups = "drop") %>%
  relocate(specific, .after = de)
cat("\nsignature genes (de = padj < 0.05 & |lfc| >= 0.25; specific = adult AUROC gate; ortholog = OMA 1:1 to chicken):\n")
print(as.data.frame(summ), row.names = FALSE)
cat("\nMAFB / EMX2 in the DE tables (excluded from every signature):\n")
print(as.data.frame(bind_rows(de) %>% filter(gene %in% TRANSGENE_TFS) %>% select(axis, gene, lfc, padj) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

write_csv(sig %>% mutate(across(c(baseMean, lfc, lfc_se, stat, adult_auroc), ~ round(.x, 4)), across(c(pvalue, padj), ~ signif(.x, 3))) %>%
            select(axis, direction, gene, chicken, lfc, padj, adult_auroc, specific, has_ortholog, excluded_tf, baseMean),
          file.path(OE_TABLES, "reference_signature_genes.csv"))
write_csv(summ, file.path(OE_TABLES, "reference_signature_summary.csv"))
iwalk(de, ~ write_csv(.x, file.path(OE_OUT, sprintf("reference_de_%s.csv.gz", .y))))
cat("\nwrote", file.path(OE_TABLES, "reference_signature_genes.csv"), "\n")
