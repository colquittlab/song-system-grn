## One clustering of the MAFB-domain accessible elements on four features, each over the lobe's full
## extent (element_profiles.py): whole-tree and oscine-lineage (Anc222) mean phyloP (accessibility-
## weighted), and mean Glut-RA and GABA-4-1 fragment signal (log1p). Each column is z-scored; lobes
## are clustered by Ward's method on Euclidean distance and drawn with their dendrogram, unsplit.
## A strip marks k groups (default 4), numbered in dendrogram order; the structure is weak (Gaussian-
## mixture BIC favours one cluster; Ward silhouettes ~0.22), so the groups are bins along gradients.
## Rows are labelled with the lobe's scored extent (1-based) and class (RA, G4, both).
##
##   Rscript plot_element_cluster.R [k=4] [depth=T] [noacc]
## thr splits rows by the rating-calibrated threshold classes (element_threshold_classes.py; implies
## noacc), ordered within each by Ward on the two features, z-scored on all lobes so colours match
## across depth filters; outputs fig_element_cluster_thr_noacc[_depth<T>].
## noacc clusters on (and shows) the two phyloP means only, naming groups from their centroids
## (deep / oscine when above average on that track, else weak); outputs carry _noacc.
## depth=T keeps only lobes with no low-depth patch: every 50-bp window of the lobe aligned in >= T of
## oscine and of non-oscine genomes (osc_depth_min50, nonosc_depth_min50; element_profiles.py);
## outputs carry _depth<T>.
## Writes results/mafb_domain/fig_element_cluster_k<k>.{pdf,png} and element_cluster_k<k>.tsv.
suppressMessages({library(ComplexHeatmap); library(circlize); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 4L
dmin <- sub("^depth=", "", grep("^depth=", args, value = TRUE)); dmin <- if (length(dmin)) as.numeric(dmin) else NA
noacc <- "noacc" %in% args                         # cluster on the two phyloP means only
thr <- "thr" %in% args                               # split rows by threshold classes (element_threshold_classes.py)
if (thr) noacc <- TRUE
tag <- paste0(if (thr) "_thr" else "", if (noacc) "_noacc" else "", if (is.na(dmin)) "" else sprintf("_depth%g", dmin))
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
feat <- c(whole_ext_mean = "whole-tree mean phyloP", osc_ext_mean = "oscine mean phyloP",
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
THR <- c("deep+oscine", "deep", "oscine", "weak")
if (thr) S[, group := thr_class]
if (!thr) print(data.table(group = nm, n = as.vector(table(cl)), round(cz, 2)))
if (!thr) fwrite(S[, .(element, cls, group, chrom, ext_start, ext_end, n_lobes, lobe, whole_ext_mean, osc_ext_mean,
             acc_ra_ext, acc_g4_ext)], file.path(res, sprintf("element_cluster_k%d%s.tsv", k, tag)), sep = "\t")

lbl <- sprintf("%s:%s-%s%s  %s", S$chrom, format(S$ext_start + 1, big.mark = ",", trim = TRUE),
               format(S$ext_end, big.mark = ",", trim = TRUE),
               ifelse(S$n_lobes > 1, sprintf(" (lobe %d/%d)", S$lobe, S$n_lobes), ""), S$cls)
gcol <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#7a3ca0", "#eda100", "#c8161a", "#5f5f5f", "#9a9a9a")[seq_len(k)], nm)
if (thr) nm <- THR
if (noacc) gcol <- setNames(ifelse(grepl("^deep\\+", nm), "#7a3ca0", ifelse(grepl("^deep", nm), "#c8161a",
                            ifelse(grepl("^oscine", nm), "#1baf7a", "#9a9a9a"))), nm)
ra <- rowAnnotation(group = S$group, class = S$cls,
                    col = list(group = gcol, class = c(RA = FIG_PAL[[2]], G4 = FIG_PAL[[1]], both = FIG_PAL[[3]])),
                    annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "bottom", annotation_name_rot = 45, simple_anno_size = unit(2, "mm"),
                    annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
ht <- Heatmap(X, name = "z", col = colorRamp2(c(-2, 0, 2), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH)),
              cluster_rows = if (thr) TRUE else hc, clustering_method_rows = "ward.D2", cluster_row_slices = FALSE,
              row_split = if (thr) factor(S$group, levels = THR) else NULL, row_title_gp = gpar(fontsize = 6), row_title_rot = 0,
              row_gap = unit(1, "mm"), row_dend_width = unit(18, "mm"), cluster_columns = FALSE,
              right_annotation = ra, row_labels = lbl, row_names_gp = gpar(fontsize = 5, fontfamily = "mono"),
              column_names_gp = gpar(fontsize = 6), column_names_rot = 45, width = unit(ncol(X) * 6, "mm"),
              column_title = if (thr) NULL else paste0(sprintf("Ward, z-scored (k = %d)", k), if (is.na(dmin)) "" else sprintf(", depth >= %g", dmin)), column_title_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
h <- 1.4 + nrow(S) * 2 / 25.4 + if (thr) 0.5 else 0
for (dev in c("pdf", "png")) {
  f <- file.path(res, if (thr) sprintf("fig_element_cluster%s.%s", tag, dev) else sprintf("fig_element_cluster_k%d%s.%s", k, tag, dev))
  if (dev == "pdf") cairo_pdf(f, width = 5.2, height = h, family = FIG_FONT) else png(f, width = 5.2, height = h, units = "in", res = 600, family = FIG_FONT)
  if (thr) {                                         # figure-level title: the two columns are too narrow for it
    draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 8, 2), "mm"))
    upViewport(0)
    grid.text(sprintf("thresholds: whole-tree mean >= %+.2f, oscine mean >= %+.2f%s", S$tcut_whole[1], S$tcut_oscine[1],
                      if (is.na(dmin)) "" else sprintf("; depth >= %g", dmin)),
              x = unit(0.5, "npc"), y = unit(1, "npc") - unit(4, "mm"), gp = gpar(fontsize = 6))
  } else draw(ht, merge_legend = TRUE)
  dev.off()
}
message("wrote ", file.path(res, if (thr) sprintf("fig_element_cluster%s", tag) else sprintf("fig_element_cluster_k%d%s", k, tag)), ".pdf/.png")
