## Constraint classes of RA and G4 accessible elements at MAFB (DHX35-TOP1 domain).
##
## Reads results/mafb_domain/element_summaries.tsv and element_profiles_{whole,oscine,access}.tsv.gz
## (element_profiles.py): SnapATAC2 peaks overlapping an own-contrast ArchR DAR, split into lobes,
## centred on accessibility summits. Lobes are clustered (Ward, Euclidean) on z-scored summaries
## over the summit +/-100 bp -- whole-tree and oscine-lineage (Anc222) phyloP: % > 1, mean, block %
## (continuity), % < -1 (acceleration) -- and drawn as
## Classes: RA, G4, and both (an RA and a G4 peak overlapping each other, merged).
##   1. a summary heatmap split by cluster, annotated with class, Log2FC and multi-part status
##   2. profile heatmaps in the same row order: whole-tree | oscine phyloP | accessibility difference
## Cluster composition by class is tested with Fisher's exact test.
##
##   Rscript plot_element_heatmaps.R [k=4] [byclass] [rule] [labels]
## byclass: rows split first by class (RA, G4, both), then by the same (shared) clusters within each.
## rule: rows split by the rule-based classes of element_rule_classes.py (deep+oscine, deep, oscine,
## weak; full-extent scores, thresholds calibrated on the visual ratings) in place of the clusters.
## labels: profiles only, one 2-mm row per lobe labelled with its scored extent (1-based, browser
## coordinates), lobe number and rule class; written as ..._labels.{pdf,png}.
## mb: rows split by the per-tree mean + block % classes (element_rule_classes.py, column mb) and
## ordered within each by Ward clustering on those four full-extent features.
## Writes results/mafb_domain/fig_element_heatmap_{summary,profiles}_k<k>[_byclass].{pdf,png} and
## element_clusters_k<k>.tsv.
suppressMessages({library(ComplexHeatmap); library(circlize); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))
## DAR classes as in the comparative_genomics figures (Glut-RA red, GABA-4 blue #4393C3; shared = their comparator
## grey); phyloP on a purple (conserved) - white - brown (accelerated) scale, clear of both class colours
CLS_COL <- c(RA = "#e31a1c", G4 = "#4393C3", both = "#9AA3AF")
PHY_DIV <- c("#b35806", "#f7f7f7", "#542788")

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args) >= 1) as.integer(args[1]) else 4L
byclass <- "byclass" %in% args; tag <- if (byclass) "_byclass" else ""
res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
S <- fread(file.path(res, "element_rule_classes.tsv"))   # element_summaries.tsv plus rule classes
rule <- "rule" %in% args; if (rule) tag <- paste0(tag, "_rule")
labels <- "labels" %in% args                        # coordinate row labels, tall layout for browsing
mb <- "mb" %in% args; if (mb) tag <- paste0(tag, "_mb")
P <- lapply(c(whole = "whole", oscine = "oscine", access = "access"), function(t) {
  m <- as.matrix(read.delim(gzfile(file.path(res, sprintf("element_profiles_%s.tsv.gz", t))), check.names = FALSE)[, -1])
  rownames(m) <- S$element; m })

feat <- c(whole_pct1 = "whole % > 1", whole_mean = "whole mean", whole_block = "whole block %", whole_acc = "whole % < -1",
          osc_pct1 = "oscine % > 1", osc_mean = "oscine mean", osc_block = "oscine block %", osc_acc = "oscine % < -1")
if (mb) feat <- c(whole_ext_mean = "whole mean", whole_ext_block = "whole block %",
                 osc_ext_mean = "oscine mean", osc_ext_block = "oscine block %")   # full extent, acc.-weighted
Xc <- scale(as.matrix(S[, names(feat), with = FALSE]))   # row-order features (summit +/-100 bp; mb: extent)
if (rule && !mb) feat <- c(feat, whole_ext_block = "whole block % (extent)", osc_ext_pct1 = "oscine % > 1 (extent)")
X <- scale(as.matrix(S[, names(feat), with = FALSE])); colnames(X) <- feat; rownames(X) <- S$element
hc <- hclust(dist(Xc), method = "ward.D2")
cl <- cutree(hc, k = k)
## number clusters by mean whole-tree constraint so the order reads high to low
ord <- order(-tapply(S$whole_pct1, cl, mean)); cl <- match(cl, ord)
S[, cluster := cl]
RULES <- c("deep+oscine", "deep", "oscine", "weak")
grp <- if (rule) match(S$rule, RULES) else S$cluster; glab <- if (rule) RULES else paste0("C", seq_len(k))
if (mb) { grp <- match(S$mb, RULES); glab <- RULES }   # per-tree mean + block % clusters (element_rule_classes.py)

## composition and per-cluster summaries
tab <- table(cluster = S$cluster, class = S$cls)
ft <- if (all(dim(tab) > 1)) fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value else NA
summ <- S[, c(list(n = .N, RA = sum(cls == "RA"), G4 = sum(cls == "G4"), both = sum(cls == "both")),
              lapply(.SD, function(v) round(mean(v), 2))), by = cluster, .SDcols = c(names(feat), "log2fc")][order(cluster)]
print(summ); message(sprintf("class x cluster Fisher p = %.3g (k = %d, %d lobes)", ft, k, nrow(S)))
if (!mb) fwrite(S[, .(element, cls, cluster, rule, mb, chrom, summit, peak_start, peak_end, ext_start, ext_end, n_lobes, lobe, forced,
             log2fc, log2fc_ra, log2fc_g4, fdr, whole_pct1, whole_mean, whole_block, whole_acc, osc_pct1, osc_mean, osc_block,
             osc_acc, whole_ext_pct1, whole_ext_block, whole_ext_mean, osc_ext_pct1, osc_ext_block, osc_ext_mean, thr_whole, thr_oscine)],
       file.path(res, sprintf("element_clusters_k%d.tsv", k)), sep = "\t")

## annotations
cls_col <- CLS_COL
ra <- rowAnnotation(class = S$cls, `Log2FC` = S$log2fc, `multi-part` = ifelse(S$n_lobes > 1, "yes", "no"),
                    col = list(class = cls_col, `Log2FC` = colorRamp2(c(0, 2, 4), c("#f0efec", "#9a9a9a", "#3a3a3a")),
                               `multi-part` = c(yes = "#5f5f5f", no = "#e6e5e1")),
                    annotation_name_gp = gpar(fontsize = 5), annotation_name_side = "top", simple_anno_size = unit(2, "mm"),
                    annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
split <- factor(glab[grp], levels = glab)
## rows: by cluster, then in dendrogram order within each cluster
rord <- hc$order[order(grp[hc$order])]
if (byclass) {
  ## class first (RA, G4, both), then the shared clusters; label slices with their sizes
  lv <- as.vector(t(outer(c("RA", "G4", "both"), glab, paste)))
  key <- paste(S$cls, glab[grp]); n <- table(key)
  lv <- lv[lv %in% names(n)]
  split <- factor(key, levels = lv, labels = sprintf("%s (%d)", lv, n[lv]))
  rord <- hc$order[order(match(S$cls[hc$order], c("RA", "G4", "both")), grp[hc$order])]
}
ht_opt$TITLE_PADDING <- unit(1, "mm")

h1 <- Heatmap(X, name = "z", col = colorRamp2(c(-2, 0, 2), PHY_DIV),
              cluster_rows = FALSE, row_order = rord, row_split = split, cluster_columns = FALSE,
              show_row_names = FALSE, column_names_gp = gpar(fontsize = 5), row_title_gp = gpar(fontsize = 6), row_title_rot = 0, row_gap = unit(if (byclass) 0.8 else 1, "mm"),
              left_annotation = ra, column_title = if (mb) "full extent, accessibility-weighted (clustering features)" else if (rule) "constraint: summit ±100 bp | full extent (rule features)" else "constraint at the accessibility summit (±100 bp)",
              column_title_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
save_ht <- function(ht, stem, w, h) {
  for (dev in c("pdf", "png")) {
    f <- file.path(res, paste0(stem, ".", dev))
    if (dev == "pdf") cairo_pdf(f, width = w, height = h, family = FIG_FONT) else png(f, width = w, height = h, units = "in", res = 600, family = FIG_FONT)
    draw(ht, merge_legend = TRUE); dev.off()
  }
  message("wrote ", file.path(res, stem), ".pdf/.png")
}
save_ht(h1, if (rule || mb) paste0("fig_element_heatmap_summary", tag) else sprintf("fig_element_heatmap_summary_k%d%s", k, tag), 3.8, if (byclass) 6.2 else 4.6)

## profiles in the same order
pos <- as.integer(colnames(P$whole)); lab <- ifelse(pos %% 250 == 0, as.character(pos), "")
mk <- function(m, nm, lim, title, cols, left = NULL) Heatmap(pmax(pmin(m, lim[2]), lim[1]), name = nm,
  col = if (length(cols) == 3) colorRamp2(c(lim[1], 0, lim[2]), cols) else colorRamp2(c(0, lim[2]), cols), cluster_rows = FALSE, row_order = rord, row_split = split,
  cluster_columns = FALSE, show_row_names = FALSE, column_labels = lab, column_names_gp = gpar(fontsize = 5),
  column_names_rot = 0, column_title = title, column_title_gp = gpar(fontsize = 6), row_title_gp = gpar(fontsize = 6), row_title_rot = 0,
  use_raster = TRUE, raster_quality = 3, width = unit(28, "mm"), left_annotation = left,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5),
                              at = if (length(cols) == 3) c(lim[1], 0, lim[2]) else c(0, lim[2]),
                              labels = round(if (length(cols) == 3) c(lim[1], 0, lim[2]) else c(0, lim[2]), 1)))
div <- PHY_DIV
alim <- as.numeric(quantile(abs(P$access), 0.98))
h2 <- mk(P$whole, "whole phyloP", c(-6, 6), "whole tree", div, left = ra) +
      mk(P$oscine, "oscine phyloP", c(-2, 2), "oscine lineage (Anc222)", div) +
      mk(P$access, "access. diff", c(0, alim), "accessibility (up − ref)", c(FIG_DIV_MID, "#2b2b2b"))
if (labels) {
  lbl <- sprintf("%s:%s-%s%s  %s  %s", S$chrom, format(S$ext_start + 1, big.mark = ",", trim = TRUE),
                 format(S$ext_end, big.mark = ",", trim = TRUE), ifelse(S$n_lobes > 1, sprintf(" (lobe %d/%d)", S$lobe, S$n_lobes), ""),
                 S$cls, if (mb) S$mb else S$rule)
  h2 <- h2 + rowAnnotation(coords = anno_text(lbl, gp = gpar(fontsize = 5, fontfamily = "mono")))
  save_ht(h2, paste0("fig_element_heatmap_profiles", tag, "_labels"), 9.0, 1.4 + nrow(S) * 2 / 25.4)
  quit(save = "no")
}
save_ht(h2, if (rule || mb) paste0("fig_element_heatmap_profiles", tag) else sprintf("fig_element_heatmap_profiles_k%d%s", k, tag), 6.4, if (byclass) 6.2 else 4.6)
