## One clustering of the MAFB-domain accessible elements on four features, each over the lobe's full
## extent (element_profiles.py): whole-tree and oscine-lineage (Anc222) mean phyloP (accessibility-
## weighted), and mean Glut-RA and GABA-4-1 fragment signal (log1p). Each column is z-scored; lobes
## are clustered by Ward's method on Euclidean distance and drawn with their dendrogram, unsplit.
## A strip marks k groups (default 4), numbered in dendrogram order; the structure is weak (Gaussian-
## mixture BIC favours one cluster; Ward silhouettes ~0.22), so the groups are bins along gradients.
## Rows are labelled with the lobe's scored extent (1-based) and class (RA, G4, both).
##
##   Rscript plot_element_cluster.R [k=4] [depth=T] [noacc]
## horiz (implies thr, nolabels) lays the elements out as columns, 5 in wide; outputs carry _horiz.
## nolabels drops the row labels and the group column, widens the class column, sets the class titles
## vertically and narrows the dendrogram (8 mm) and phyloP cells (4 mm) (outputs carry _nolabels).
## thr splits rows by the rating-calibrated threshold classes (element_threshold_classes.py; implies
## noacc), ordered within each by Ward on the two features, z-scored on all lobes so colours match
## across depth filters; outputs fig_element_cluster_thr_noacc[_depth<T>].
## noacc clusters on the two phyloP means only (z-scored) and shows their actual values, each on its
## own diverging colour scale (avian / whole tree +/-5, oscine +/-1), naming groups from their centroids
## (deep / oscine when above average on that track, else weak); outputs carry _noacc.
## depth=T keeps only lobes with no low-depth patch: every 50-bp window of the lobe aligned in >= T of
## oscine and of non-oscine genomes (osc_depth_min50, nonosc_depth_min50; element_profiles.py);
## outputs carry _depth<T>.
## Writes results/mafb_domain/fig_element_cluster_k<k>.{pdf,png} and element_cluster_k<k>.tsv.
suppressMessages({library(ComplexHeatmap); library(circlize); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))
## DAR classes as in the comparative_genomics figures (Glut-RA red, GABA-4 blue #4393C3; shared = their comparator
## grey); phyloP on a purple (conserved) - white - brown (accelerated) scale, clear of both class colours
CLS_COL <- c(RA = "#e31a1c", G4 = "#4393C3", both = "#9AA3AF")
PHY_DIV <- c("#b35806", "#f7f7f7", "#542788")

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 4L
dmin <- sub("^depth=", "", grep("^depth=", args, value = TRUE)); dmin <- if (length(dmin)) as.numeric(dmin) else NA
noacc <- "noacc" %in% args                         # cluster on the two phyloP means only
thr <- "thr" %in% args                               # split rows by threshold classes (element_threshold_classes.py)
nolab <- "nolabels" %in% args                        # no row labels or group column; wider class column
horiz <- "horiz" %in% args                           # elements as columns, 5 in wide (thr, noacc)
fig6 <- "fig6" %in% args                             # Fig. 6 panel D as laid out in the figure (implies thr, nolabels)
if (horiz) { thr <- TRUE; nolab <- TRUE }
if (fig6) { thr <- TRUE; nolab <- TRUE }
if (thr) noacc <- TRUE
tag <- paste0(if (thr) "_thr" else "", if (noacc) "_noacc" else "", if (is.na(dmin)) "" else sprintf("_depth%g", dmin),
              if (nolab) "_nolabels" else "", if (horiz) "_horiz" else "", if (fig6) "_fig6" else "")
res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
S <- fread(file.path(res, "element_rule_classes.tsv"))
S0 <- copy(S)                                        # full set, for fixed z-scaling in thr mode
if (thr) S <- merge(S, fread(file.path(res, "element_threshold_classes.tsv"))[, .(element, thr_class, tcut_whole = thr_whole, tcut_oscine = thr_oscine)], by = "element")
if (!is.na(dmin)) {                                  # confidently aligned lobes only
  n0 <- nrow(S); S <- S[osc_depth_min50 >= dmin & nonosc_depth_min50 >= dmin]
  message(sprintf("depth >= %g across the lobe (every 50-bp window, both clades) keeps %d of %d lobes (%s)", dmin, nrow(S), n0,
                  paste(names(table(S$cls)), table(S$cls), collapse = ", ")))
}
S[, `:=`(lacc_ra = log1p(acc_ra_ext), lacc_g4 = log1p(acc_g4_ext))]
feat <- c(whole_ext_mean = "avian mean phyloP", osc_ext_mean = "oscine mean phyloP",
          lacc_ra = "RA accessibility", lacc_g4 = "G4 accessibility")
if (noacc) feat <- feat[1:2]
X <- scale(as.matrix(S[, names(feat), with = FALSE])); colnames(X) <- feat
if (thr) {                                           # z-score on all lobes, so colours match across depth filters
  X0 <- as.matrix(S0[, names(feat), with = FALSE])
  X <- sweep(sweep(as.matrix(S[, names(feat), with = FALSE]), 2, colMeans(X0)), 2, apply(X0, 2, sd), "/"); colnames(X) <- feat
}
hc <- hclust(dist(X), method = "ward.D2")
cl <- cutree(hc, k = k)
## number groups in dendrogram order (top to bottom)
cl <- match(cl, unique(cl[hc$order])); nm <- paste0("C", seq_len(k)); S[, group := nm[cl]]
cz <- apply(X, 2, function(v) tapply(v, cl, mean))
if (noacc) {                                         # name from centroids: above-average whole tree / oscine
  nm <- ifelse(cz[, 1] > 0 & cz[, 2] > 0, "deep+oscine", ifelse(cz[, 1] > 0, "deep", ifelse(cz[, 2] > 0, "oscine", "weak")))
  nm <- make.unique(nm, sep = " "); S[, group := nm[cl]]
}
THR <- c("deep+oscine", "strong avian", "moderate avian", "oscine", "weak")
if (thr) S[, group := thr_class]
if (!thr) print(data.table(group = nm, n = as.vector(table(cl)), round(cz, 2)))
if (!thr) fwrite(S[, .(element, cls, group, chrom, ext_start, ext_end, n_lobes, lobe, whole_ext_mean, osc_ext_mean,
             acc_ra_ext, acc_g4_ext)], file.path(res, sprintf("element_cluster_k%d%s.tsv", k, tag)), sep = "\t")

lbl <- sprintf("%s:%s-%s%s  %s", S$chrom, format(S$ext_start + 1, big.mark = ",", trim = TRUE),
               format(S$ext_end, big.mark = ",", trim = TRUE),
               ifelse(S$n_lobes > 1, sprintf(" (lobe %d/%d)", S$lobe, S$n_lobes), ""), S$cls)
gcol <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#7a3ca0", "#eda100", "#c8161a", "#5f5f5f", "#9a9a9a")[seq_len(k)], nm)
if (thr) nm <- THR
if (noacc) gcol <- setNames(ifelse(grepl("^deep\\+", nm), "#d4a017", ifelse(grepl("^strong", nm), "#3f007d", ifelse(grepl("^(deep|moderate)", nm), "#8c6bb1",
                            ifelse(grepl("^oscine", nm), "#1baf7a", "#bdbdbd")))), nm)   # as plot_element_class_bars.R
## display labels from here on (Fig. 6): DAR classes RA / PVALB-1 / shared, and the class "avian+oscine"
DISP <- c(RA = "RA", G4 = "PVALB-1", both = "shared")
S[, cls := factor(DISP[cls], levels = DISP)]
CLS_COL <- setNames(CLS_COL[names(DISP)], DISP)
if (thr) { S[group == "deep+oscine", group := "avian+oscine"]; THR[THR == "deep+oscine"] <- "avian+oscine" }
CCOL <- CLS_COL
ra <- if (fig6) rowAnnotation(class = S$cls, col = list(class = CCOL), show_annotation_name = FALSE, show_legend = FALSE,
                    gp = gpar(col = "white", lwd = 0.6), simple_anno_size = unit(5.5, "mm")) else
      if (nolab) rowAnnotation(class = S$cls, col = list(class = CCOL),
                    annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "bottom", annotation_name_rot = 45, simple_anno_size = unit(4, "mm"),
                    annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5))) else
      rowAnnotation(group = S$group, class = S$cls,
                    col = list(group = gcol, class = CCOL),
                    annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "bottom", annotation_name_rot = 45, simple_anno_size = unit(2, "mm"),
                    annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
if (noacc) {
  ## the actual means, each track on its own diverging scale; rows still
  ## ordered by Ward on the z-scored pair (clustering_distance_rows looks the z-scores up by row name)
  Xr <- as.matrix(S[, names(feat), with = FALSE]); rownames(Xr) <- rownames(X) <- S$element
  zdist <- function(m) dist(X[rownames(m), , drop = FALSE])
  LIM <- c(5, 1)                                     # avian +/-5, oscine +/-1
  one <- function(j, nm, extra = list()) do.call(Heatmap, c(list(Xr[, j, drop = FALSE], name = nm,
    col = colorRamp2(c(-LIM[j], 0, LIM[j]), PHY_DIV), cluster_columns = FALSE,
    column_labels = if (fig6) c("avian", "oscine")[j] else feat[j], column_names_gp = gpar(fontsize = if (fig6) 7 else 6),
    column_names_rot = if (fig6) 0 else 45, column_names_centered = fig6, width = unit(if (fig6) 7 else if (nolab) 4 else 6, "mm"),
    show_heatmap_legend = !fig6,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5), at = c(-LIM[j], 0, LIM[j]))), extra))
  ht <- one(1, "avian mean phyloP", list(
          cluster_rows = if (thr) TRUE else hc, clustering_distance_rows = zdist, clustering_method_rows = "ward.D2",
          cluster_row_slices = FALSE, row_split = if (thr) factor(S$group, levels = THR) else NULL,
          row_title_gp = gpar(fontsize = if (fig6) 7 else 6), row_title_rot = if (nolab) 90 else 0, row_gap = unit(1, "mm"),
          row_dend_width = unit(if (nolab) 8 else 18, "mm"),
          show_row_names = FALSE,
          column_title = if (thr) NULL else paste0(sprintf("Ward, z-scored (k = %d)", k), if (is.na(dmin)) "" else sprintf(", depth >= %g", dmin)),
          column_title_gp = gpar(fontsize = 6))) +
        one(2, "oscine mean phyloP", list(right_annotation = ra, row_labels = lbl, show_row_names = !nolab,
          row_names_gp = gpar(fontsize = 5, fontfamily = "mono")))
} else {
ht <- Heatmap(X, name = "z", col = colorRamp2(c(-2, 0, 2), PHY_DIV),
              cluster_rows = if (thr) TRUE else hc, clustering_method_rows = "ward.D2", cluster_row_slices = FALSE,
              row_split = if (thr) factor(S$group, levels = THR) else NULL, row_title_gp = gpar(fontsize = 6), row_title_rot = 0,
              row_gap = unit(1, "mm"), row_dend_width = unit(18, "mm"), cluster_columns = FALSE,
              right_annotation = ra, row_labels = lbl, row_names_gp = gpar(fontsize = 5, fontfamily = "mono"),
              column_names_gp = gpar(fontsize = 6), column_names_rot = 45, width = unit(ncol(X) * 6, "mm"),
              column_title = if (thr) NULL else paste0(sprintf("Ward, z-scored (k = %d)", k), if (is.na(dmin)) "" else sprintf(", depth >= %g", dmin)), column_title_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
}
h <- 1.4 + nrow(S) * 2 / 25.4 + if (thr) 0.5 else 0
W <- if (nolab) 3.2 else if (noacc) 6.0 else 5.2
if (nolab) h <- min(h, 5)                             # no row labels: rows shrink to keep the figure <= 5 in tall
if (fig6) W <- 2.6                                   # legends drawn compactly beside the columns (fig6_legends)
if (horiz) {
  ## elements as columns split by threshold class (Ward within each on the z-scored pair), the two phyloP
  ## means as rows, DAR class beneath; legends along the bottom so the figure is 5 in wide
  Xt <- t(Xr)
  cdist <- function(m) {                             # ComplexHeatmap passes columns as rows here
    ids <- if (!is.null(rownames(m)) && all(rownames(m) %in% rownames(X))) rownames(m) else colnames(m)
    dist(X[ids, , drop = FALSE])
  }
  lp <- function(at) list(direction = "horizontal", title_position = "topcenter", legend_width = unit(16, "mm"),
                          title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5), at = at)
  ba <- HeatmapAnnotation(class = S$cls, col = list(class = CCOL), simple_anno_size = unit(3, "mm"),
                          annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "left",
                          annotation_legend_param = list(class = list(title = "DAR class", nrow = 1, title_position = "topcenter",
                                                                      title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5))))
  mkh <- function(j, nm, extra = list()) do.call(Heatmap, c(list(Xt[j, , drop = FALSE], name = nm,
    col = colorRamp2(c(-LIM[j], 0, LIM[j]), PHY_DIV), cluster_rows = FALSE, row_labels = c("avian", "oscine")[j], row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, height = unit(5, "mm"),
    heatmap_legend_param = lp(c(-LIM[j], 0, LIM[j]))), extra))
  ht <- mkh(1, "avian\nmean phyloP", list(cluster_columns = TRUE, clustering_distance_columns = cdist,
          clustering_method_columns = "ward.D2", cluster_column_slices = FALSE, column_split = factor(S$group, levels = THR),
          column_title_gp = gpar(fontsize = 6), column_gap = unit(1, "mm"), column_dend_height = unit(4, "mm"))) %v%
        mkh(2, "oscine\nmean phyloP", list(bottom_annotation = ba))
  ht_opt$TITLE_PADDING <- unit(1, "mm")               # group labels close to the dendrogram
  h <- 2.1; W <- 5
}
for (dev in c("pdf", "png")) {
  f <- file.path(res, if (thr) sprintf("fig_element_cluster%s.%s", tag, dev) else sprintf("fig_element_cluster_k%d%s.%s", k, tag, dev))
  if (dev == "pdf") cairo_pdf(f, width = W, height = h, family = FIG_FONT) else png(f, width = W, height = h, units = "in", res = 600, family = FIG_FONT)
  if (fig6) {
    ## no title; "phyloP" over the two phyloP columns and "DARs" over the class column; colour bars and the DAR key
    ## drawn beside the columns: limit above and below each bar, its title rotated alongside
    draw(ht, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, padding = unit(c(2, 2, 6, 17), "mm"))
    decorate_heatmap_body("avian mean phyloP", slice = 1, {
      grid.text("phyloP", x = unit(1, "npc") + unit(1, "mm"), y = unit(1, "npc") + unit(1.2, "mm"), vjust = 0, gp = gpar(fontsize = 7))
    })
    decorate_annotation("class", slice = 1, {
      grid.text("DARs", y = unit(1, "npc") + unit(1.2, "mm"), vjust = 0, gp = gpar(fontsize = 7))
    })
    upViewport(0)
    bar <- function(ytop, lim, title) {               # ytop: npc of the page
      x0 <- unit(1, "npc") - unit(15, "mm"); bw <- unit(3.8, "mm"); bh <- unit(20, "mm")
      ramp <- colorRamp2(c(-lim, 0, lim), PHY_DIV)(seq(lim, -lim, length.out = 100))
      grid.raster(matrix(ramp, ncol = 1), x = x0, y = unit(ytop, "npc") - unit(3, "mm"), width = bw, height = bh,
                  just = c("left", "top"), interpolate = TRUE)
      grid.text(sprintf("%g", lim), x = x0 + bw * 0.5, y = unit(ytop, "npc") - unit(2.6, "mm"), vjust = 0, gp = gpar(fontsize = 6))
      grid.text(sprintf("%g", -lim), x = x0 + bw * 0.5, y = unit(ytop, "npc") - unit(3.4, "mm") - bh, vjust = 1, gp = gpar(fontsize = 6))
      grid.text(title, x = x0 + bw + unit(1.6, "mm"), y = unit(ytop, "npc") - unit(3, "mm") - bh * 0.5, rot = -90,
                gp = gpar(fontsize = 6))
    }
    bar(0.86, 5, "avian mean phyloP"); bar(0.57, 1, "oscine mean phyloP")
    for (i in seq_along(CLS_COL))
      grid.text(names(CLS_COL)[i], x = unit(1, "npc") - unit(16, "mm"), y = unit(0.33, "npc") - unit((i - 1) * 3.4, "mm"),
                just = "left", gp = gpar(fontsize = 8, fontface = "bold", col = CLS_COL[i]))
  } else if (thr) {                                  # figure-level title: the two columns are too narrow for it
    if (horiz) draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", legend_gap = unit(8, "mm"),
                    padding = unit(c(2, 2, 8, 2), "mm"))
    else draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 8, 2), "mm"))
    upViewport(0)
    grid.text(sprintf(if (horiz) "thresholds: avian mean >= %+.2f (strong >= 3), oscine mean >= %+.2f%s"
                      else if (nolab) "avian mean >= %+.2f (strong >= 3)\noscine mean >= %+.2f%s"
                      else "thresholds: avian mean >= %+.2f (strong >= 3), oscine mean >= %+.2f%s", S$tcut_whole[1], S$tcut_oscine[1],
                      if (is.na(dmin)) "" else sprintf("; depth >= %g", dmin)),
              x = unit(0.5, "npc"), y = unit(1, "npc") - unit(4, "mm"), gp = gpar(fontsize = 6))
  } else draw(ht, merge_legend = TRUE)
  dev.off()
}
message("wrote ", file.path(res, if (thr) sprintf("fig_element_cluster%s", tag) else sprintf("fig_element_cluster_k%d%s", k, tag)), ".pdf/.png")
