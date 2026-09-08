## ---------------------------------------------------------------------------
## Composite / per-method heatmaps for the two full-suite hybrid-label
## cross-species composites:
##   finch x chicken  (gg_adult_hybrid,  Zaremba et al. adult chicken)
##   finch x mouse    (yao_adult_hybrid, Yao et al. 2023 ABC Atlas)
##
## R port of finch-integration-toolkit/plot_rank_heatmap.py (ComplexHeatmap in
## place of seaborn.clustermap), so these figures are built with the same
## toolchain (renv + config/figure_theme.R) as the rest of the repo. Design
## decisions carried over unchanged:
##   * MATRIX COLOUR: a [0,1] magnitude gets a single-hue sequential ramp
##     (Oranges, fixed 0-1 so panels are comparable across comparisons); a
##     signed correlation (GSI) gets a blue-white-orange diverging ramp on [-1,1].
##   * CLASS STRIPS: coarse class of every finch cluster (row) and reference
##     label (column), six display classes; identity is also in the legend.
##   * COLUMN FILTERING: a reference has 50-500 labels; columns are the union of
##     each finch cluster's top-k matches (+ each row's peak-agreement column).
##     The number dropped is stated in the subtitle, never silently.
##   * DOTS: N dots in a cell = N methods independently calling that pair a
##     reciprocal top-N match (Zaremba convention); drawn for N >= dots_min.
##   * GEOMETRY: square cells whose pitch derives from the tick-label size, so
##     labels pack tightly; --scale shrinks the whole figure proportionally.
##
## Inputs are the tracked CSVs written by assemble_{gg,yao}_adult_hybrid.py into
## composite_scoring/results/<tag>/ plus the reference-label annotation tables in
## composite_scoring/annotations/. Outputs (PDF + PNG, gitignored) overwrite the
## files of the same name in each results dir.
##
## Usage:
##   Rscript snrna/integration/plot_composite_heatmaps_hybrid.R [tag ...] [name ...]
##   (no args = every figure for both tags)
## ---------------------------------------------------------------------------

suppressMessages({
  library(here)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
source(here::here("config/figure_theme.R"))
source(here::here("snrna/integration/R/composite_heatmap_utils.R"))

## DIMNAME_PADDING is specifically the gap between the heatmap body and the
## row/column NAME text (not the class-strip annotations -- ROW_ANNO_PADDING /
## COLUMN_ANNO_PADDING govern those, left at the default 1mm). Cut as close to
## 0 as ComplexHeatmap allows without the text touching the cells.
ht_opt$DIMNAME_PADDING <- unit(0, "mm")

RES <- here::here("snrna/integration/composite_scoring/results")
ANN <- here::here("snrna/integration/composite_scoring/annotations")
ANNOT <- c(gg_adult_hybrid  = file.path(ANN, "gg_adult_label_annotation.csv"),
           yao_adult_hybrid = file.path(ANN, "yao_label_annotation.csv"))

LABEL_PT <- 6.5          # tick-label size; the square cell pitch derives from it
MIN_FONT_PT <- 6         # no rendered text may go below this, max_width_in yields first
## Cell pitch (row/column height/width) as a multiple of the label font size in
## points -- this, not DIMNAME_PADDING, is what sets the spacing between one row
## label and the next (row names are one line, not rotated, so adjacent labels
## collide if pitch drops below roughly the font's own line height). 1.05 is
## close to the floor for that: single-line text at fontsize pt has a line
## height of about 1.0-1.15x pt with this font, so much less than 1.0 starts
## clipping ascenders/descenders between rows.
PITCH_MULT <- 1.05
ORANGES <- c("#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c",
             "#f16913", "#d94801", "#a63603", "#7f2704")

## --- Renderer -------------------------------------------------------------------
## `scale` below is only the STARTING guess for the cell/font/legend scale
## factor; the actual factor used is solved for so the drawn figure (row labels
## included) is never wider than `max_width_in` -- see the auto-fit loop at the
## bottom of this function.
plot_rank_heatmap <- function(matrix_csv, out_prefix, annot_csv = NULL, top_k = 3,
                              max_cols = 90, title = "composite rank score",
                              cbar_label = "composite rank score", signed = FALSE,
                              dots_csv = NULL, dots_min = 2, transpose = FALSE,
                              scale = 1, label_pt_override = NULL, max_width_in = 7) {
  M <- read_matrix(matrix_csv)
  annot <- if (!is.null(annot_csv) && file.exists(annot_csv))
    read.csv(annot_csv, row.names = 1, check.names = FALSE) else NULL
  n_all <- ncol(M)
  dots_full <- NULL
  if (!is.null(dots_csv) && file.exists(dots_csv)) {
    dd <- read_matrix(dots_csv)
    dots_full <- matrix(0L, nrow(M), ncol(M), dimnames = dimnames(M))
    ri <- intersect(rownames(dd), rownames(M)); ci <- intersect(colnames(dd), colnames(M))
    dots_full[ri, ci] <- as.integer(round(dd[ri, ci]))
  }
  keep <- filter_columns(M, top_k, max_cols, dots_full, dots_min)
  Msub <- M[, keep, drop = FALSE]
  dropped <- n_all - length(keep)
  message(sprintf("%s: %d x %d -> kept %d columns (union of per-cluster top-%d); dropped %d",
                  basename(matrix_csv), nrow(M), n_all, length(keep), top_k, dropped))
  dots <- if (!is.null(dots_full)) dots_full[rownames(Msub), keep, drop = FALSE] else NULL

  ## classes are assigned BEFORE any transpose (finch rule for finch clusters,
  ## reference rule for reference labels), so a flip is purely presentational
  row_cls <- display_class(rownames(Msub), expected_class(rownames(Msub)))
  col_cls <- display_class(colnames(Msub), ref_class(colnames(Msub), annot))
  if (transpose) {
    Msub <- t(Msub); if (!is.null(dots)) dots <- t(dots)
    tmp <- row_cls; row_cls <- col_cls; col_cls <- tmp
  }

  ## hierarchical clustering: euclidean / average, as the Python tool
  hc_row <- hclust(dist(Msub), method = "average")
  hc_col <- hclust(dist(t(Msub)), method = "average")

  vmin <- if (signed) -1 else 0
  col_fun <- if (signed) colorRamp2(c(-1, 0, 1), c("#2a78d6", "#ffffff", "#eb6834"))
             else colorRamp2(seq(0, 1, length.out = length(ORANGES)), ORANGES)

  present <- CLASS_ORDER[CLASS_ORDER %in% union(row_cls, col_cls)]
  if ("unknown" %in% union(row_cls, col_cls)) present <- c(present, "unknown")
  cls_col <- CLASS_COLORS[present]

  scale_note <- if (signed) "-1 to 1 (diverging, signed)" else "0-1"
  ## one clause per line: the heatmap is only pitch*ncol wide, a single line clips
  subtitle <- paste(c(
    sprintf("%s — colour scale fixed %s for cross-panel comparability; observed range here %.2f–%.2f",
            gsub("\n", " ", cbar_label), scale_note, min(Msub), max(Msub)),
    sprintf("%d of %d reference labels shown (union of per-cluster top-%d; %d omitted)",
            length(keep), n_all, top_k, dropped),
    if (!is.null(dots)) sprintf("dots = methods agreeing this pair is a reciprocal top-N match (Zaremba); shown for ≥%d", dots_min)),
    collapse = "\n")

  ## Every size below is expressed in terms of `sc`, so a single scalar controls
  ## the whole figure's footprint -- what the auto-fit loop searches over to hit
  ## max_width_in without a second, independent knob to keep in sync. Every
  ## fontsize is floored at MIN_FONT_PT regardless of how far `sc` shrinks, so
  ## max_width_in is the soft constraint here, not the font floor -- a wide
  ## matrix can end up wider than max_width_in rather than illegible.
  build <- function(sc) {
    label_pt <- max(LABEL_PT * sc, MIN_FONT_PT)
    render_pt <- max(if (is.null(label_pt_override)) label_pt else label_pt_override, MIN_FONT_PT)
    pitch <- label_pt / 72 * PITCH_MULT             # inches per row / column (square)
    strip_in <- 0.09 * sc
    gp_lab <- gpar(fontsize = render_pt, fontfamily = FIG_FONT)

    left_anno <- rowAnnotation(class = row_cls, col = list(class = cls_col),
                               show_annotation_name = FALSE, show_legend = FALSE,
                               simple_anno_size = unit(strip_in, "in"))
    top_anno <- HeatmapAnnotation(class = col_cls, col = list(class = cls_col),
                                  show_annotation_name = FALSE, show_legend = FALSE,
                                  simple_anno_size = unit(strip_in, "in"))

    dot_d <- pitch * 0.20                           # dot diameter, inches
    cell_fun <- NULL
    if (!is.null(dots)) {
      cell_fun <- function(j, i, x, y, w, h, fill) {
        n <- dots[i, j]
        if (n >= dots_min) {
          off <- dot_offsets(n)
          grid.points(x + unit(off[, 1] * pitch, "in"), y + unit(off[, 2] * pitch, "in"),
                      pch = 16, size = unit(dot_d, "in"), gp = gpar(col = "black"))
        }
      }
    }

    ht <- Heatmap(
      Msub, name = "score", col = col_fun,
      cluster_rows = hc_row, cluster_columns = hc_col,
      show_row_dend = FALSE, show_column_dend = FALSE,
      width = unit(pitch * ncol(Msub), "in"), height = unit(pitch * nrow(Msub), "in"),
      row_names_gp = gp_lab, column_names_gp = gp_lab, column_names_rot = 90,
      row_names_side = "right", left_annotation = left_anno, top_annotation = top_anno,
      cell_fun = cell_fun, use_raster = FALSE, border = FALSE,
      column_title = paste0(title, "\n", subtitle),
      column_title_gp = gpar(fontsize = max(7.5 * sc, MIN_FONT_PT), fontfamily = FIG_FONT),
      heatmap_legend_param = list(
        title = cbar_label, at = c(vmin, 1), labels = c(vmin, 1),
        title_gp = gpar(fontsize = max(6.5 * sc, MIN_FONT_PT), fontfamily = FIG_FONT),
        labels_gp = gpar(fontsize = max(6 * sc, MIN_FONT_PT), fontfamily = FIG_FONT),
        legend_height = unit(1.1 * sc, "in"), grid_width = unit(0.08 * sc, "in"))
    )
    class_lgd <- Legend(labels = present, legend_gp = gpar(fill = cls_col),
                        title = "coarse class (colour strips)", ncol = min(length(present), 4L), by_row = TRUE,
                        title_gp = gpar(fontsize = max(7.5 * sc, MIN_FONT_PT), fontfamily = FIG_FONT),
                        labels_gp = gpar(fontsize = max(7 * sc, MIN_FONT_PT), fontfamily = FIG_FONT),
                        grid_height = unit(0.1 * sc, "in"), grid_width = unit(0.18 * sc, "in"))
    list(ht = ht, class_lgd = class_lgd)
  }

  draw_built <- function(b) draw(b$ht, annotation_legend_list = list(b$class_lgd),
                                 annotation_legend_side = "top", heatmap_legend_side = "right",
                                 merge_legend = FALSE, padding = unit(c(2, 2, 2, 2), "mm"))

  ## Measure the drawn object on a throwaway cairo device (cairo, not pdf(): the
  ## latter has no Arial in its PostScript font database and warns per string).
  measure <- function(sc) {
    tmp <- tempfile(fileext = ".pdf")
    cairo_pdf(tmp, width = 30, height = 30, family = FIG_FONT)
    d <- draw_built(build(sc))
    w <- convertWidth(ComplexHeatmap:::width(d), "in", valueOnly = TRUE)
    h <- convertHeight(ComplexHeatmap:::height(d), "in", valueOnly = TRUE)
    dev.off(); unlink(tmp)
    c(w, h)
  }

  ## Auto-fit: shrink sc (from the job's own starting scale) until the drawn
  ## width -- row labels, strips, dendrogram and legend included, since that's
  ## everything measure() returns -- is at most max_width_in. Every term in
  ## build() scales with sc except the 4mm fixed padding, so one proportional
  ## step gets close and a couple more converge; sc only ever shrinks here, never
  ## grows past the job's own starting point. Stops at sc_min_font, the point
  ## where build()'s own MIN_FONT_PT floor kicks in -- past that, shrinking sc
  ## further does not shrink the (floored) label font or the pitch derived from
  ## it, so it would only spin iterations without changing the drawn width.
  sc_min_font <- MIN_FONT_PT / LABEL_PT
  sc <- scale
  wh <- measure(sc)
  iter <- 0L
  while (wh[1] > max_width_in + 0.02 && iter < 6L && sc > sc_min_font + 1e-6) {
    sc <- max(sc * max_width_in / wh[1], sc_min_font)
    wh <- measure(sc)
    iter <- iter + 1L
  }
  w <- wh[1]; h <- wh[2]
  if (w > max_width_in + 0.02)
    message(sprintf("  NOTE: %s stays %.1fin wide (> max_width_in %.1f) -- MIN_FONT_PT %d reached first",
                    basename(out_prefix), w, max_width_in, MIN_FONT_PT))

  fig_check_font()
  b <- build(sc)
  cairo_pdf(paste0(out_prefix, ".pdf"), width = w, height = h, family = FIG_FONT)
  draw_built(b); dev.off()
  png(paste0(out_prefix, ".png"), width = w, height = h, units = "in", res = 220,
      type = "cairo", family = FIG_FONT, bg = "white")
  draw_built(b); dev.off()
  message(sprintf("  wrote %s.pdf / .png  (%.1fx%.1f in, scale %.2f)", out_prefix, w, h, sc))
  invisible(c(w, h))
}

## --- Figure catalogue ---------------------------------------------------------------
METHOD_TITLES <- list(gsi = list("GSI (correlation)", TRUE),
                      samap = list("SAMap alignment score", FALSE),
                      cca = list("Seurat CCA transfer score", FALSE),
                      saturn = list("SATURN transfer score", FALSE))
GSI_C2021 <- c(all = "GSI, Colquitt 2021 method (DEG markers, all genes)",
               nontf = "GSI, Colquitt 2021 method (variable genes, non-TF)",
               tf = "GSI, Colquitt 2021 method (variable genes, TF-only)")

jobs <- function(tag) {
  d <- file.path(RES, tag)
  conf <- list(matrix_csv = file.path(d, "composite_confidence_matrix.csv"),
               cbar_label = "mapping confidence",
               dots_csv = file.path(d, "composite_agreement_count_matrix.csv"))
  J <- list(rank_score_clustermap = list(matrix_csv = file.path(d, "composite_rank_score_matrix.csv"),
                                         title = "Composite rank-aggregate score"))
  if (tag == "gg_adult_hybrid") {
    J$confidence_clustermap <- c(conf, title = "Composite mapping confidence")
    for (k in names(GSI_C2021))
      J[[sprintf("method_gsi_colquitt2021_%s_heatmap", k)]] <- list(
        matrix_csv = file.path(RES, sprintf("gsi_corr_gg_adult_hybrid_colquitt2021_%s.csv", k)),
        title = GSI_C2021[[k]], signed = TRUE)
  } else {
    t <- list(title = "Composite mapping confidence (finch x Yao mouse)")
    J$confidence_clustermap <- c(conf, t)
    J$confidence_clustermap_transposed <- c(conf, t, transpose = TRUE, scale = 0.8)
    J$confidence_clustermap_topk1_anchored <- c(conf, t, top_k = 1)
    J$confidence_clustermap_topk1_anchored_transposed <- c(conf, t, top_k = 1, transpose = TRUE)
    for (pt in c(5, 6))
      J[[sprintf("confidence_clustermap_topk1_anchored_transposed_pt%d", pt)]] <-
        c(conf, t, top_k = 1, transpose = TRUE, scale = 0.6, label_pt_override = pt)
  }
  for (m in names(METHOD_TITLES))
    J[[sprintf("method_%s_heatmap", m)]] <- list(
      matrix_csv = file.path(d, sprintf("method_%s_matrix.csv", m)),
      title = METHOD_TITLES[[m]][[1]], signed = METHOD_TITLES[[m]][[2]])
  J
}

main <- function(argv) {
  tags <- intersect(argv, names(ANNOT)); if (!length(tags)) tags <- names(ANNOT)
  only <- setdiff(argv, names(ANNOT))
  for (tag in tags) {
    J <- jobs(tag)
    for (nm in names(J)) {
      if (length(only) && !nm %in% only) next
      message(sprintf("\n>>> %s/%s", tag, nm))
      do.call(plot_rank_heatmap, c(J[[nm]], list(out_prefix = file.path(RES, tag, nm),
                                                 annot_csv = ANNOT[[tag]])))
    }
  }
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
