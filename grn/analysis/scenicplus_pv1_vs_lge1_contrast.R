#!/usr/bin/env Rscript
# RNA contrast PVALB-1 vs GABA-LGE-1 for every gene, as an alternative comparator to the MAFB-low MGE type LAMP5 (scenicplus_mafb_pvalb_vs_mafb_low_mge.R).
# Same method as there: pseudobulk (sample_short x assignment x cluster, >= MIN_CELLS cells), DESeq2 on counts, design ~ library + group when estimable.
# The original LAMP5-based table (scenicplus_all_genes_onoff_contrasts_config37.csv.gz) is not touched. Output: scenicplus_all_genes_pv1_vs_lge1_config37.csv.gz
# with gene, lfc_PV1_vs_LGE1, padj_PV1_vs_LGE1 (positive = higher in PVALB-1).
#
#   Rscript scenicplus_pv1_vs_lge1_contrast.R
suppressMessages({library(Seurat); library(qs2); library(DESeq2); library(Matrix); library(tidyverse); library(here)})
set.seed(2026)
MAIN <- "config37"
OBJ <- "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
HERE <- here::here("grn/analysis")
PV1 <- "GABA-MGE-PVALB-1"; LGE1 <- "GABA-LGE-1"; MIN_CELLS <- 10
obj <- qs_read(OBJ, nthreads = 8)
meta <- obj@meta.data; meta$cluster <- meta$cluster_hybrid
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
rm(obj); invisible(gc())
grp <- paste(meta$sample_short, meta$assignment, meta$cluster, sep = "|")
f <- factor(grp)
pb <- cnt %*% sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
colnames(pb) <- levels(f)
parts <- do.call(rbind, strsplit(levels(f), "|", fixed = TRUE))
info <- tibble(id = levels(f), lib = parts[, 1], cluster = parts[, 3], n = as.integer(table(f)))
s <- info %>% filter(cluster %in% c(PV1, LGE1), n >= MIN_CELLS) %>% mutate(g = ifelse(cluster == PV1, "fg", "bg"))
cd <- data.frame(row.names = s$id, lib = factor(s$lib), grp = factor(s$g, levels = c("bg", "fg")))
nfg <- sum(cd$grp == "fg"); nbg <- sum(cd$grp == "bg"); stopifnot(nfg >= 2, nbg >= 2)
mm <- model.matrix(~ lib + grp, cd)
design <- if (qr(mm)$rank == ncol(mm)) ~ lib + grp else ~ grp
dds <- DESeqDataSetFromMatrix(round(as.matrix(pb[, s$id])), cd, design)
dds <- DESeq(dds[rowSums(counts(dds) >= 5) >= 2, ], quiet = TRUE)
r <- results(dds, contrast = c("grp", "fg", "bg"))
cat(sprintf("PVALB-1 vs LGE-1: pseudobulk samples %d vs %d (cells: %d vs %d) | design %s | genes tested %d\n", nfg, nbg,
            sum(meta$cluster == PV1), sum(meta$cluster == LGE1), format(design), nrow(r)))
out <- tibble(gene = rownames(r), lfc_PV1_vs_LGE1 = round(r$log2FoldChange, 3), padj_PV1_vs_LGE1 = signif(r$padj, 3))
write_csv(out, file.path(HERE, paste0("scenicplus_all_genes_pv1_vs_lge1_", MAIN, ".csv.gz")))
for (g in c("MAFB", "PVALB", "KCNC1", "ERBB4")) cat(sprintf("  %-6s log2FC %5.2f padj %s\n", g, out$lfc_PV1_vs_LGE1[out$gene == g], signif(out$padj_PV1_vs_LGE1[out$gene == g], 2)))
