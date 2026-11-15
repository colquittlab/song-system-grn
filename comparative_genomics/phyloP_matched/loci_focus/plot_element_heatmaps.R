## Constraint classes of RA and G4 accessible elements at MAFB (DHX35-TOP1 domain).
##
## Reads results/mafb_domain/element_summaries.tsv and element_profiles_{whole,oscine,access}.tsv.gz
## (element_profiles.py): SnapATAC2 peaks overlapping an own-contrast ArchR DAR, split into lobes,
## centred on accessibility summits. Lobes are clustered (Ward, Euclidean) on z-scored summaries
## over the summit +/-100 bp -- whole-tree and oscine-lineage (Anc222) phyloP: % > 1, mean, block %
## (continuity), % < -1 (acceleration) -- and drawn as
##   1. a summary heatmap split by cluster, annotated with class, Log2FC and multi-part status
##   2. profile heatmaps in the same row order: whole-tree | oscine phyloP | accessibility difference
## Cluster composition by class is tested with Fisher's exact test.
##
##   Rscript plot_element_heatmaps.R [k=4] [byclass]
## byclass: rows split first by class (RA, G4), then by the same (shared) clusters within each.
## Writes results/mafb_domain/fig_element_heatmap_{summary,profiles}_k<k>[_byclass].{pdf,png} and
## element_clusters_k<k>.tsv.
suppressMessages({library(ComplexHeatmap); library(circlize); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args) >= 1) as.integer(args[1]) else 4L
byclass <- "byclass" %in% args; tag <- if (byclass) "_byclass" else ""
res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
S <- fread(file.path(res, "element_summaries.tsv"))
P <- lapply(c(whole = "whole", oscine = "oscine", access = "access"), function(t) {
  m <- as.matrix(read.delim(gzfile(file.path(res, sprintf("element_profiles_%s.tsv.gz", t))), check.names = FALSE)[, -1])
  rownames(m) <- S$element; m })

feat <- c(whole_pct1 = "whole % > 1", whole_mean = "whole mean", whole_block = "whole block %", whole_acc = "whole % < -1",
          osc_pct1 = "oscine % > 1", osc_mean = "oscine mean", osc_block = "oscine block %", osc_acc = "oscine % < -1")
X <- scale(as.matrix(S[, names(feat), with = FALSE])); colnames(X) <- feat; rownames(X) <- S$element
hc <- hclust(dist(X), method = "ward.D2")
cl <- cutree(hc, k = k)
## number clusters by mean whole-tree constraint so the order reads high to low
ord <- order(-tapply(S$whole_pct1, cl, mean)); cl <- match(cl, ord)
S[, cluster := cl]

## composition and per-cluster summaries
tab <- table(cluster = S$cluster, class = S$cls)
ft <- if (all(dim(tab) > 1)) fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value else NA
summ <- S[, c(list(n = .N, RA = sum(cls == "RA"), G4 = sum(cls == "G4")),
              lapply(.SD, function(v) round(mean(v), 2))), by = cluster, .SDcols = c(names(feat), "log2fc")][order(cluster)]
print(summ); message(sprintf("class x cluster Fisher p = %.3g (k = %d, %d lobes)", ft, k, nrow(S)))
fwrite(S[, .(element, cls, cluster, chrom, summit, peak_start, peak_end, n_lobes, lobe, forced, log2fc, fdr,
             whole_pct1, whole_mean, whole_block, whole_acc, osc_pct1, osc_mean, osc_block, osc_acc)],
       file.path(res, sprintf("element_clusters_k%d.tsv", k)), sep = "\t")

## annotations
cls_col <- c(RA = FIG_PAL[[2]], G4 = FIG_PAL[[1]])
ra <- rowAnnotation(class = S$cls, `Log2FC` = S$log2fc, `multi-part` = ifelse(S$n_lobes > 1, "yes", "no"),
                    col = list(class = cls_col, `Log2FC` = colorRamp2(c(0, 2, 4), c("#f0efec", "#9a9a9a", "#3a3a3a")),
                               `multi-part` = c(yes = "#5f5f5f", no = "#e6e5e1")),
                    annotation_name_gp = gpar(fontsize = 5), annotation_name_side = "top", simple_anno_size = unit(2, "mm"),
                    annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
split <- factor(paste0("C", S$cluster), levels = paste0("C", seq_len(k)))
## rows: by cluster, then in dendrogram order within each cluster
rord <- hc$order[order(S$cluster[hc$order])]
if (byclass) {
  ## class first (RA, G4), then the shared clusters; label slices with their sizes
  lv <- as.vector(t(outer(c("RA", "G4"), paste0("C", seq_len(k)), paste)))
  key <- paste(S$cls, paste0("C", S$cluster)); n <- table(key)
  lv <- lv[lv %in% names(n)]
  split <- factor(key, levels = lv, labels = sprintf("%s (%d)", lv, n[lv]))
  rord <- hc$order[order(match(S$cls[hc$order], c("RA", "G4")), S$cluster[hc$order])]
}
ht_opt$TITLE_PADDING <- unit(1, "mm")

h1 <- Heatmap(X, name = "z", col = colorRamp2(c(-2, 0, 2), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH)),
              cluster_rows = FALSE, row_order = rord, row_split = split, cluster_columns = FALSE,
              show_row_names = FALSE, column_names_gp = gpar(fontsize = 5), row_title_gp = gpar(fontsize = 6), row_title_rot = 0, row_gap = unit(if (byclass) 0.8 else 1, "mm"),
              left_annotation = ra, column_title = "constraint at the accessibility summit (±100 bp)",
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
save_ht(h1, sprintf("fig_element_heatmap_summary_k%d%s", k, tag), 3.8, if (byclass) 5.6 else 4.6)

## profiles in the same order
pos <- as.integer(colnames(P$whole)); lab <- ifelse(pos %% 250 == 0, as.character(pos), "")
mk <- function(m, nm, lim, title, cols, left = NULL) Heatmap(pmax(pmin(m, lim[2]), lim[1]), name = nm,
  col = if (length(cols) == 3) colorRamp2(c(lim[1], 0, lim[2]), cols) else colorRamp2(c(0, lim[2]), cols), cluster_rows = FALSE, row_order = rord, row_split = split,
  cluster_columns = FALSE, show_row_names = FALSE, column_labels = lab, column_names_gp = gpar(fontsize = 5),
  column_names_rot = 0, column_title = title, column_title_gp = gpar(fontsize = 6), row_title_gp = gpar(fontsize = 6), row_title_rot = 0,
  use_raster = TRUE, raster_quality = 3, width = unit(28, "mm"), left_annotation = left,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
div <- c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH)
alim <- as.numeric(quantile(abs(P$access), 0.98))
h2 <- mk(P$whole, "whole phyloP", c(-3, 3), "whole tree", div, left = ra) +
      mk(P$oscine, "oscine phyloP", c(-2, 2), "oscine lineage (Anc222)", div) +
      mk(P$access, "access. diff", c(0, alim), "accessibility (up − ref)", c(FIG_DIV_MID, "#2b2b2b"))
save_ht(h2, sprintf("fig_element_heatmap_profiles_k%d%s", k, tag), 6.4, if (byclass) 5.6 else 4.6)
