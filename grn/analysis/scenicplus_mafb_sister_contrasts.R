#!/usr/bin/env Rscript
# MAFB targets and PAIRWISE (sister) differences between related cell types.
#
# scenicplus_mafb_shared_program.R scores each gene by its specificity across ALL neuron types. That can hide a gene that is clearly
# higher in RA than in its sister C1H-1 but also expressed in many other neurons. Here each focal type (RA, PVALB-1, PVALB-2, SST-1)
# is compared with its nearest sister, with pseudobulk DESeq2 (replicate = library x bird x cell type, >= 10 cells; design ~ library +
# type, so dissection/library effects are modeled when the two types share libraries). The sister is the nearest other neuron cluster
# by correlation of mean expression over the 3,000 most variable genes (printed, so "C1H-1 is RA's sister" is checked, not assumed).
#
# For the MAFB +/+ targets of the main config (config37):
#   - how many are up in each focal type vs its sister (padj < 0.05, log2FC > 0.25), against random genes matched on expression level and
#     breadth (so a regulon of broadly expressed genes is not credited);
#   - whether the genes up in RA vs C1H-1 are also up in an interneuron type vs ITS sister (count in both, and the correlation of
#     the two log2FCs), against the same matched null;
#   - what the cross-cell-type specificity score of the earlier analysis said about the RA-vs-sister genes.
# Caveat built into the design: PVALB-1 and PVALB-2 come largely from different dissections (hvc/nc vs arco), so their sister contrast
# leans on the few libraries containing both; RA vs C1H-1 is paired within the two RA libraries.
#
#   Rscript scenicplus_mafb_sister_contrasts.R
suppressMessages({library(Seurat); library(qs2); library(DESeq2); library(Matrix); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
source(here::here("multiome/archr/hybrid_labels.R"))   # hybrid_ct_order
set.seed(2026)
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
OBJ <- "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
HERE <- here::here("grn/analysis")
FOCAL <- c(RA = "Glut-CACNA1H-RA", PV1 = "GABA-MGE-PVALB-1", PV2 = "GABA-MGE-PVALB-2", SST1 = "GABA-MGE-SST-1")
LFC <- 0.25; PADJ <- 0.05; MIN_CELLS <- 10; NPERM <- 2000
LABEL <- c("PVALB", "KCNC1", "MAFB", "PCSK5", "PIK3R6", "ABHD12", "FAM184B", "RUNX2", "PIK3R5", "ADAMTS18", "MTNR1B", "ST6GALNAC5",
           "PLXNA4", "TENM4", "KIAA1217", "CNTNAP2", "VIPR2", "ERBB4", "ARX", "RARB", "COL4A1", "SGCZ", "LOC110473696")

obj <- qs_read(OBJ, nthreads = 8)
obj$cluster <- obj$cluster_hybrid
meta <- obj@meta.data
n_cl <- table(meta$cluster)
avg <- as.matrix(log1p(AverageExpression(obj, assays = "SCT", layer = "data", group.by = "cluster")[[1]]))
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
rm(obj); invisible(gc())

## 1. sisters from the data
neu <- intersect(hybrid_ct_order, names(n_cl)[n_cl >= 50]); neu <- grep("^(Glut|GABA)-", neu, value = TRUE)
vg <- names(sort(apply(avg[, neu], 1, var), decreasing = TRUE))[1:3000]
cc <- cor(avg[vg, neu])
sister <- sapply(FOCAL, function(f) { o <- sort(cc[f, setdiff(neu, f)], decreasing = TRUE); names(o)[1] })
cat("nearest neuron clusters (correlation over 3,000 variable genes):\n")
for (f in FOCAL) { o <- sort(cc[f, setdiff(neu, f)], decreasing = TRUE); cat(sprintf("  %-18s %s\n", f, paste0(names(o)[1:3], " (", round(o[1:3], 2), ")", collapse = ", "))) }

## 2. pseudobulk
grp <- paste(meta$sample_short, meta$assignment, meta$cluster, sep = "|")
f <- factor(grp)
M <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
pb <- cnt %*% M
colnames(pb) <- levels(f)
parts <- do.call(rbind, strsplit(levels(f), "|", fixed = TRUE))
info <- tibble(id = levels(f), lib = parts[, 1], asg = parts[, 2], cluster = parts[, 3], n = as.integer(table(f)))

de <- function(fg, bg) {
  s <- info %>% filter(cluster %in% c(fg, bg), n >= MIN_CELLS)
  cd <- data.frame(row.names = s$id, lib = factor(s$lib), grp = factor(s$cluster, levels = c(bg, fg)))
  nfg <- sum(cd$grp == fg); nbg <- sum(cd$grp == bg)
  stopifnot(nfg >= 2, nbg >= 2)
  mm <- model.matrix(~ lib + grp, cd)
  design <- if (qr(mm)$rank == ncol(mm)) ~ lib + grp else ~ grp
  mat <- round(as.matrix(pb[, s$id]))
  dds <- DESeqDataSetFromMatrix(mat, cd, design)
  dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, contrast = c("grp", fg, bg))
  cat(sprintf("  %-18s vs %-18s pseudobulk samples %d vs %d | design %s | genes tested %d\n", fg, bg, nfg, nbg, format(design), nrow(r)))
  tibble(gene = rownames(r), lfc = r$log2FoldChange, padj = r$padj)
}
cat("\npseudobulk DESeq2 (replicate = library x bird x cell type, >= ", MIN_CELLS, " cells):\n", sep = "")
res <- imap(FOCAL, function(fg, nm) de(fg, sister[[nm]]))
ra_user <- de(FOCAL[["RA"]], "Glut-CACNA1H-1")   # the contrast the analysis was built around, whatever the data-derived sister is

## 3. MAFB +/+ targets and the matched null
e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name))
pool <- Reduce(intersect, c(lapply(res, function(d) d$gene[!is.na(d$padj)]), list(ra_user$gene[!is.na(ra_user$padj)], rownames(avg))))
mu <- rowMeans(avg[pool, neu]); sdv <- apply(avg[pool, neu], 1, sd)
strata <- setNames(paste(cut(mu, quantile(mu, 0:10 / 10), include.lowest = TRUE, labels = FALSE),
                         cut(sdv, quantile(sdv, 0:4 / 4), include.lowest = TRUE, labels = FALSE), sep = "_"), pool)
tg <- intersect(unique(e$Gene), pool)
cat("\n", MAIN, " MAFB +/+ targets: ", length(unique(e$Gene)), " | testable in all contrasts: ", length(tg), " | matching pool: ", length(pool), " genes\n", sep = "")
tab <- function(nm, d) d %>% filter(gene %in% pool) %>% select(gene, lfc, padj) %>% rename_with(~ paste0(.x, "_", nm), c(lfc, padj))
W <- reduce(c(imap(res, ~ tab(.y, .x)), list(tab("RA_vs_C1H1", ra_user))), full_join, by = "gene")
up <- function(nm) !is.na(W[[paste0("padj_", nm)]]) & W[[paste0("padj_", nm)]] < PADJ & W[[paste0("lfc_", nm)]] > LFC
Wm <- W %>% column_to_rownames("gene")
upv <- sapply(c(names(FOCAL), "RA_vs_C1H1"), function(nm) setNames(up(nm), W$gene))
lfcv <- sapply(c(names(FOCAL), "RA_vs_C1H1"), function(nm) setNames(W[[paste0("lfc_", nm)]], W$gene))
write_csv(W %>% mutate(across(starts_with("lfc"), ~ round(.x, 3)), across(starts_with("padj"), ~ signif(.x, 3))) %>% filter(gene %in% tg),
          file.path(HERE, paste0("scenicplus_mafb_sister_contrasts_", MAIN, ".csv")))

by_s <- function(excl) split(setdiff(pool, excl), strata[setdiff(pool, excl)])
draw <- function(need, bs) unlist(lapply(names(need), function(k) { p <- bs[[k]]; if (!length(p)) character() else sample(p, need[[k]], replace = length(p) < need[[k]]) }))
need <- table(strata[tg]); bs <- by_s(tg)
nulls <- replicate(NPERM, draw(need, bs), simplify = FALSE)
ratio <- function(obs, nullv) c(obs = obs, null = mean(nullv), ratio = obs / mean(nullv), p = (1 + sum(nullv >= obs)) / (length(nullv) + 1))
rows <- list()
for (nm in c(names(FOCAL), "RA_vs_C1H1")) {
  r <- ratio(sum(upv[tg, nm]), sapply(nulls, function(s) sum(upv[s, nm])))
  rows[[length(rows) + 1]] <- tibble(test = "targets up vs sister", pair = nm, obs = r["obs"], null = r["null"], ratio = r["ratio"], p = r["p"])
}
for (nm in c("RA", "RA_vs_C1H1")) for (pp in c("PV1", "PV2", "SST1")) {
  both <- function(g) sum(upv[g, nm] & upv[g, pp])
  r <- ratio(both(tg), sapply(nulls, both))
  rows[[length(rows) + 1]] <- tibble(test = "targets up in both", pair = paste(nm, pp, sep = " + "), obs = r["obs"], null = r["null"], ratio = r["ratio"], p = r["p"])
  cr <- function(g) suppressWarnings(cor(lfcv[g, nm], lfcv[g, pp], use = "complete.obs", method = "spearman"))
  o <- cr(tg); nv <- sapply(nulls, cr)
  rows[[length(rows) + 1]] <- tibble(test = "log2FC Spearman", pair = paste(nm, pp, sep = " + "), obs = o, null = mean(nv), ratio = (o - mean(nv)) / sd(nv), p = (1 + sum(nv >= o)) / (length(nv) + 1))
}
S <- bind_rows(rows)
write_csv(S %>% mutate(across(where(is.numeric), ~ round(.x, 3))), file.path(HERE, paste0("scenicplus_mafb_sister_enrichment_", MAIN, ".csv")))
cat("\n(ratio = observed / matched-null mean; for 'log2FC Spearman' the column is the z-score vs the null)\n")
print(as.data.frame(S %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

## 4. what did the cross-cell-type score say about the genes that ARE up in RA vs its sister?
gs <- read_csv(file.path(OUT, "mafb_gene_specificity_RA_PVALB2.csv"), show_col_types = FALSE)
ra_up <- intersect(tg[upv[tg, "RA_vs_C1H1"]], gs$gene)
zz <- gs$z_RA[match(ra_up, gs$gene)]
cat(sprintf("\nMAFB targets up in RA vs C1H-1: %d (of %d testable); their cross-neuron specificity z in RA: median %.2f; z > 1 for %d (%.0f%%), z <= 0.5 for %d (%.0f%%)\n",
            sum(upv[tg, "RA_vs_C1H1"]), length(tg), median(zz), sum(zz > 1), 100 * mean(zz > 1), sum(zz <= 0.5), 100 * mean(zz <= 0.5)))
for (pp in c("PV1", "PV2", "SST1")) {
  b <- tg[upv[tg, "RA_vs_C1H1"] & upv[tg, pp]]
  cat(sprintf("  up in RA vs C1H-1 AND in %s vs its sister (%s): %d genes: %s\n", FOCAL[[pp]], sister[[pp]], length(b), paste(head(b, 25), collapse = ", ")))
}

## 5. heatmap: sister-contrast log2FC of the targets that are up in at least one
cols <- c(RA_vs_C1H1 = "RA vs C1H-1", RA = paste0("RA vs ", sub("^Glut-", "", sister[["RA"]])), PV1 = paste0("PVALB-1 vs ", sub("^GABA-MGE-", "", sister[["PV1"]])),
          PV2 = paste0("PVALB-2 vs ", sub("^GABA-MGE-", "", sister[["PV2"]])), SST1 = paste0("SST-1 vs ", sub("^GABA-MGE-", "", sister[["SST1"]])))
if (sister[["RA"]] == "Glut-CACNA1H-1") cols <- cols[names(cols) != "RA"]       # same contrast, show once
if (sister[["PV1"]] == FOCAL[["PV2"]] && sister[["PV2"]] == FOCAL[["PV1"]]) cols <- cols[names(cols) != "PV2"]   # mutual sisters: mirror image
sel <- tg[rowSums(upv[tg, names(cols), drop = FALSE]) > 0]
L <- lfcv[sel, names(cols), drop = FALSE]; L[is.na(L)] <- 0
ord <- rownames(L)[hclust(dist(L), method = "ward.D2")$order]
ZC <- 2
long <- as_tibble(lfcv[sel, names(cols), drop = FALSE], rownames = "gene") %>%
  pivot_longer(-gene, names_to = "contrast", values_to = "lfc") %>%
  left_join(as_tibble(upv[sel, names(cols), drop = FALSE], rownames = "gene") %>% pivot_longer(-gene, names_to = "contrast", values_to = "sig"), by = c("gene", "contrast")) %>%
  mutate(gene = factor(gene, levels = rev(ord)), contrast = factor(cols[contrast], levels = cols), lfcc = pmax(pmin(lfc, ZC), -ZC))
lab <- intersect(LABEL, sel)
nx <- length(cols); lab_df <- tibble(gene = factor(lab, levels = rev(ord)), x = nx + 0.8)
p <- ggplot(long, aes(contrast, gene, fill = lfcc)) +
  geom_raster() +
  ggrepel::geom_text_repel(data = lab_df, aes(x = x, y = gene, label = gene), inherit.aes = FALSE, hjust = 0, direction = "y", xlim = c(nx + 0.8, NA),
                           size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY, segment.size = 0.2, segment.color = FIG_INK_MUTED,
                           min.segment.length = 0, box.padding = 0.12, point.padding = 0, max.overlaps = Inf, seed = 1) +
  scale_fill_gradient2(low = FIG_DIV_LOW, mid = FIG_DIV_MID, high = FIG_DIV_HIGH, midpoint = 0, limits = c(-ZC, ZC), na.value = FIG_NEUTRAL_FILL,
                       name = "log2 fold change", breaks = c(-2, 0, 2), labels = c("≤−2", "0", "≥2")) +
  scale_x_discrete(expand = expansion(add = c(0, 4.5))) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FIG_PT_AXIS_TEXT), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(), legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.9, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
fig_save(p, file.path(OUT, MAIN, "mafb_target_sister_contrasts_heatmap"), width = 3, height = 5)
cat("\nheatmap rows (targets up in >= 1 sister contrast):", length(sel), "\n")
