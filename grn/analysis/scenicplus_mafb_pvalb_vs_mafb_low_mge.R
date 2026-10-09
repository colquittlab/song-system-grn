#!/usr/bin/env Rscript
# A MAFB on/off contrast for the interneuron side: PVALB-1 (and PVALB-2) against the OTHER MGE types that do not express MAFB.
#
# In scenicplus_mafb_sister_contrasts.R the nearest sister of PVALB-1 is PVALB-2, and MAFB is not differential between them, so that
# contrast cannot test whether the RA program is reused. RA vs C1H-1 is a MAFB on/off contrast (MAFB +2.1 log2FC in RA). The
# matching contrast here is PVALB vs the MGE types where MAFB is low.
#
# "MAFB-low" is derived, not assumed: every MGE cluster with >= 50 cells is compared with PVALB-1 by pseudobulk DESeq2 and the table of
# MAFB expression (mean log1p, % of cells with a count) and MAFB log2FC is printed. A cluster counts as MAFB-low if MAFB is lower than in
# PVALB-1 by log2FC > 0.5 (padj < 0.05; contrast = PVALB-1 vs the other cluster, so positive means MAFB is higher in PVALB-1). The contrast pools those clusters (design ~ library + group when estimable).
# Then, for the config37 MAFB +/+ targets: targets up in PVALB vs MAFB-low MGE, overlap with targets up in RA vs C1H-1 and the
# correlation of the two log2FCs, each against random genes matched on expression and breadth.
#
#   Rscript scenicplus_mafb_pvalb_vs_mafb_low_mge.R
suppressMessages({library(Seurat); library(qs2); library(DESeq2); library(Matrix); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
set.seed(2026)
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
OBJ <- "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
HERE <- here::here("grn/analysis")
RA <- "Glut-CACNA1H-RA"; C1H1 <- "Glut-CACNA1H-1"
PV1 <- "GABA-MGE-PVALB-1"; PV2 <- "GABA-MGE-PVALB-2"
LFC <- 0.25; PADJ <- 0.05; MIN_CELLS <- 10; NPERM <- 2000
LABEL <- c("PVALB", "KCNC1", "MAFB", "MAF", "PCSK5", "PIK3R6", "ABHD12", "FAM184B", "RUNX2", "PIK3R5", "ADAMTS18", "MTNR1B", "ST6GALNAC5",
           "PLXNA4", "TENM4", "KIAA1217", "CNTNAP2", "VIPR2", "ERBB4", "ARX", "GRIK1", "ELFN1", "MGLL", "PDZD2", "COL19A1", "AFDN")

obj <- qs_read(OBJ, nthreads = 8)
obj$cluster <- obj$cluster_hybrid
meta <- obj@meta.data
n_cl <- table(meta$cluster)
avg <- as.matrix(log1p(AverageExpression(obj, assays = "SCT", layer = "data", group.by = "cluster")[[1]]))
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
rm(obj); invisible(gc())
mge <- sort(grep("^GABA-MGE-", names(n_cl)[n_cl >= 50], value = TRUE))
neu <- grep("^(Glut|GABA)-", names(n_cl)[n_cl >= 50], value = TRUE)
cat("MGE clusters with >= 50 cells:", paste(mge, collapse = ", "), "\n")

grp <- paste(meta$sample_short, meta$assignment, meta$cluster, sep = "|")
f <- factor(grp)
pb <- cnt %*% sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
colnames(pb) <- levels(f)
parts <- do.call(rbind, strsplit(levels(f), "|", fixed = TRUE))
info <- tibble(id = levels(f), lib = parts[, 1], cluster = parts[, 3], n = as.integer(table(f)))

de <- function(fg, bg, label) {
  s <- info %>% filter(cluster %in% c(fg, bg), n >= MIN_CELLS) %>% mutate(g = ifelse(cluster %in% fg, "fg", "bg"))
  cd <- data.frame(row.names = s$id, lib = factor(s$lib), grp = factor(s$g, levels = c("bg", "fg")))
  nfg <- sum(cd$grp == "fg"); nbg <- sum(cd$grp == "bg")
  stopifnot(nfg >= 2, nbg >= 2)
  mm <- model.matrix(~ lib + grp, cd)
  design <- if (qr(mm)$rank == ncol(mm)) ~ lib + grp else ~ grp
  dds <- DESeqDataSetFromMatrix(round(as.matrix(pb[, s$id])), cd, design)
  dds <- DESeq(dds[rowSums(counts(dds) >= 5) >= 2, ], quiet = TRUE)
  r <- results(dds, contrast = c("grp", "fg", "bg"))
  cat(sprintf("  %-34s pseudobulk samples %2d vs %2d | design %-10s | genes tested %d\n", label, nfg, nbg, format(design), nrow(r)))
  tibble(gene = rownames(r), lfc = r$log2FoldChange, padj = r$padj)
}

## 1. where is MAFB? expression per MGE cluster, and PVALB-1 vs each other MGE cluster
mafb_ct <- as.numeric(as(cnt["MAFB", , drop = FALSE], "dgCMatrix")[1, ])
cat("\nMAFB in the MGE clusters (and in RA / C1H-1 for reference):\n")
ex <- tibble(cluster = c(RA, C1H1, mge), n_cells = as.integer(n_cl[c(RA, C1H1, mge)]), mean_log1p = round(avg["MAFB", c(RA, C1H1, mge)], 2),
             pct_cells = round(100 * sapply(c(RA, C1H1, mge), function(cl) mean(mafb_ct[meta$cluster == cl] > 0))))
print(as.data.frame(ex), row.names = FALSE)
cat("\npseudobulk DESeq2, PVALB-1 vs each other MGE cluster:\n")
vs <- map(setdiff(mge, PV1), function(o) de(PV1, o, paste(PV1, "vs", o)))
names(vs) <- setdiff(mge, PV1)
mafb_lfc <- map_dfr(names(vs), function(o) vs[[o]] %>% filter(gene == "MAFB") %>% mutate(other = o))
print(as.data.frame(mafb_lfc %>% mutate(lfc = round(lfc, 2), padj = signif(padj, 2)) %>% select(other, lfc, padj)), row.names = FALSE)
low <- mafb_lfc %>% filter(!is.na(padj), padj < PADJ, lfc > 0.5) %>% pull(other)
cat("\nMAFB-low MGE clusters (MAFB lower than in PVALB-1: PVALB-1 vs other log2FC > 0.5, padj < 0.05):", if (length(low)) paste(low, collapse = ", ") else "none", "\n")
stopifnot(length(low) >= 1)

## 2. contrasts
cat("\ncontrasts:\n")
C <- list(RA_vs_C1H1 = de(RA, C1H1, "RA vs C1H-1"),
          PV1_vs_low = de(PV1, low, paste("PVALB-1 vs MAFB-low MGE (", paste(sub("GABA-MGE-", "", low), collapse = "+"), ")", sep = "")))
mafb_pv2 <- de(PV2, setdiff(mge, PV2), "PVALB-2 vs all other MGE (info only)") %>% filter(gene == "MAFB")
cat("  MAFB in PVALB-2 vs all other MGE: log2FC", round(mafb_pv2$lfc, 2), "padj", signif(mafb_pv2$padj, 2), "\n")
low_for_pv2 <- setdiff(low, PV2)
C$PV2_vs_low <- de(PV2, low_for_pv2, paste("PVALB-2 vs MAFB-low MGE (", paste(sub("GABA-MGE-", "", low_for_pv2), collapse = "+"), ")", sep = ""))
for (nm in names(C)) cat(sprintf("  MAFB in %-12s log2FC %5.2f padj %s\n", nm, C[[nm]]$lfc[C[[nm]]$gene == "MAFB"], signif(C[[nm]]$padj[C[[nm]]$gene == "MAFB"], 2)))

## 3. targets, matched null
e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name))
pool <- Reduce(intersect, c(lapply(C, function(d) d$gene[!is.na(d$padj)]), list(rownames(avg))))
mu <- rowMeans(avg[pool, neu]); sdv <- apply(avg[pool, neu], 1, sd)
strata <- setNames(paste(cut(mu, quantile(mu, 0:10 / 10), include.lowest = TRUE, labels = FALSE), cut(sdv, quantile(sdv, 0:4 / 4), include.lowest = TRUE, labels = FALSE), sep = "_"), pool)
tg <- intersect(unique(e$Gene), pool)
cat("\n", MAIN, " MAFB +/+ targets: ", length(unique(e$Gene)), " | testable: ", length(tg), " | pool: ", length(pool), "\n", sep = "")
tab <- function(nm, d) d %>% filter(gene %in% pool) %>% rename_with(~ paste0(.x, "_", nm), c(lfc, padj))
W <- reduce(imap(C, ~ tab(.y, .x)), full_join, by = "gene")
upv <- sapply(names(C), function(nm) setNames(!is.na(W[[paste0("padj_", nm)]]) & W[[paste0("padj_", nm)]] < PADJ & W[[paste0("lfc_", nm)]] > LFC, W$gene))
lfcv <- sapply(names(C), function(nm) setNames(W[[paste0("lfc_", nm)]], W$gene))
write_csv(W %>% filter(gene %in% tg) %>% mutate(across(starts_with("lfc"), ~ round(.x, 3)), across(starts_with("padj"), ~ signif(.x, 3))),
          file.path(HERE, paste0("scenicplus_mafb_pvalb_vs_mafb_low_mge_", MAIN, ".csv")))
bs <- split(setdiff(pool, tg), strata[setdiff(pool, tg)]); need <- table(strata[tg])
draw <- function() unlist(lapply(names(need), function(k) { p <- bs[[k]]; if (!length(p)) character() else sample(p, need[[k]], replace = length(p) < need[[k]]) }))
nulls <- replicate(NPERM, draw(), simplify = FALSE)
rt <- function(test, pair, obs, nv, z = FALSE) tibble(test = test, pair = pair, obs = obs, null = mean(nv),
  ratio = if (z) (obs - mean(nv)) / sd(nv) else obs / mean(nv), p = (1 + sum(nv >= obs)) / (length(nv) + 1))
rows <- list()
for (nm in names(C)) rows[[length(rows) + 1]] <- rt("targets up", nm, sum(upv[tg, nm]), sapply(nulls, function(s) sum(upv[s, nm])))
for (pp in c("PV1_vs_low", "PV2_vs_low")) {
  both <- function(g) sum(upv[g, "RA_vs_C1H1"] & upv[g, pp])
  rows[[length(rows) + 1]] <- rt("targets up in both", paste("RA_vs_C1H1 +", pp), both(tg), sapply(nulls, both))
  cr <- function(g) suppressWarnings(cor(lfcv[g, "RA_vs_C1H1"], lfcv[g, pp], use = "complete.obs", method = "spearman"))
  rows[[length(rows) + 1]] <- rt("log2FC Spearman (ratio col = z)", paste("RA_vs_C1H1 +", pp), cr(tg), sapply(nulls, cr), z = TRUE)
}
## directional criterion: no padj (the interneuron contrasts, with few pseudobulk samples, are underpowered for it); log2FC above a threshold
for (thr in c(0, 0.25, 0.5, 1)) {
  for (nm in names(C)) {
    dn <- function(g) sum(lfcv[g, nm] > thr, na.rm = TRUE)
    rows[[length(rows) + 1]] <- rt(paste0("directional: log2FC > ", thr), nm, dn(tg), sapply(nulls, dn))
  }
  for (pp in c("PV1_vs_low", "PV2_vs_low")) {
    dd <- function(g) sum(lfcv[g, "RA_vs_C1H1"] > thr & lfcv[g, pp] > thr, na.rm = TRUE)
    rows[[length(rows) + 1]] <- rt(paste0("directional: both log2FC > ", thr), paste("RA_vs_C1H1 +", pp), dd(tg), sapply(nulls, dd))
  }
}
S <- bind_rows(rows)
write_csv(S %>% mutate(across(where(is.numeric), ~ round(.x, 3))), file.path(HERE, paste0("scenicplus_mafb_pvalb_vs_mafb_low_mge_enrichment_", MAIN, ".csv")))
cat("\n"); print(as.data.frame(S %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
for (pp in c("PV1_vs_low", "PV2_vs_low")) {
  b <- tg[upv[tg, "RA_vs_C1H1"] & upv[tg, pp]]
  cat(sprintf("\nup in RA vs C1H-1 AND in %s: %d genes: %s\n", pp, length(b), paste(b, collapse = ", ")))
}

## 4. heatmap: log2FC in the MAFB on/off contrasts, targets up in at least one
cols <- c(RA_vs_C1H1 = "RA vs C1H-1", PV1_vs_low = "PVALB-1 vs MAFB-low MGE", PV2_vs_low = "PVALB-2 vs MAFB-low MGE")
sel <- tg[rowSums(upv[tg, names(cols), drop = FALSE]) > 0]
L <- lfcv[sel, names(cols), drop = FALSE]; L[is.na(L)] <- 0
ord <- rownames(L)[hclust(dist(L), method = "ward.D2")$order]
ZC <- 2
long <- as_tibble(lfcv[sel, names(cols), drop = FALSE], rownames = "gene") %>% pivot_longer(-gene, names_to = "contrast", values_to = "lfc") %>%
  mutate(gene = factor(gene, levels = rev(ord)), contrast = factor(cols[contrast], levels = cols), lfcc = pmax(pmin(lfc, ZC), -ZC))
lab <- intersect(LABEL, sel); nx <- length(cols)
lab_df <- tibble(gene = factor(lab, levels = rev(ord)), x = nx + 0.8)
p <- ggplot(long, aes(contrast, gene, fill = lfcc)) + geom_raster() +
  ggrepel::geom_text_repel(data = lab_df, aes(x = x, y = gene, label = gene), inherit.aes = FALSE, hjust = 0, direction = "y", xlim = c(nx + 0.8, NA),
                           size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, segment.color = FIG_INK_MUTED,
                           min.segment.length = 0, box.padding = 0.12, point.padding = 0, max.overlaps = Inf, seed = 1) +
  scale_fill_gradient2(low = FIG_DIV_LOW, mid = FIG_DIV_MID, high = FIG_DIV_HIGH, midpoint = 0, limits = c(-ZC, ZC), na.value = FIG_NEUTRAL_FILL,
                       name = "log2 fold change", breaks = c(-2, 0, 2), labels = c("≤−2", "0", "≥2")) +
  scale_x_discrete(expand = expansion(add = c(0, 3.5))) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_fig() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FIG_PT_AXIS_TEXT), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(), legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.9, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
fig_save(p, file.path(OUT, MAIN, "mafb_target_mafb_onoff_contrasts_heatmap"), width = 3, height = 5)
cat("\nheatmap rows:", length(sel), "\n")
