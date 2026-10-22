#!/usr/bin/env Rscript
# Figure 2, panel B -- integration of song and non-song glutamatergic neurons: UMAPs before and after integration, and the reciprocal
# label-transfer score of each song type with each non-song type ("integration partner").
#
# Song (HVC, LMAN, RA types) and non-song glutamatergic neurons are treated as two datasets. Each is normalized on its own with SCTransform
# (the object's SCT assay, fitted on all cells together, already contains the song/non-song difference), with 3000 integration features
# (SelectIntegrationFeatures, PrepSCTIntegration) and a 30-PC PCA per dataset.
#   UMAPs:   left, the two datasets merged without correction (PCA on the integration features, UMAP on PCs 1-30); right, after CCA
#            integration (FindIntegrationAnchors, k.anchor 20, dims 1-30; IntegrateData; PCA and UMAP on PCs 1-30 of the integrated assay).
#            Seurat defaults otherwise (UMAP: uwot, cosine, 30 neighbours, min.dist 0.3, seed 42). Song types are coloured by nucleus,
#            non-song types grey; song-type labels at each type's median position.
#   Heatmap: label transfer in both directions (FindTransferAnchors, CCA, k.anchor 20, dims 1-30; TransferData, k.weight 50). Per direction
#            and query type a, reference type b: M[a, b] = (share of a's cells predicted as b) x (mean prediction score for b over a's cells).
#            Reciprocal score S[song, non-song] = sqrt(fwd[song, non-song] x rev[non-song, song]); high only if each type's cells vote for
#            the other. Colour scale from 0 to the panel's maximum.
#
# This is the only tracked copy of this code. It reproduces three notebooks that are not in the repository:
#   snrna/integration/song_nonsong_integration_hybrid.qmd          per-dataset SCTransform preps (cache/prep_ALL.qs2)
#   snrna/integration/song_nonsong_umap_hybrid.qmd                 -> song_nonsong_umap_hybrid/umap_all_by_type_bare.pdf      (UMAPs)
#   snrna/integration/song_nonsong_reciprocal_transfer_hybrid.qmd  -> song_nonsong_reciprocal_transfer_hybrid/reciprocal_all_cca_k20_bare.pdf
# using their "ALL" analysis (every glutamatergic type) at CCA, k.anchor 20.
#
# Input: the glutamatergic snRNA-seq object with the hybrid labels (common.R). The per-dataset prep is cached in FIG2_OUT/panel_B_prep_ALL.qs2
#        and reused if present; FIG2_B_PREP points at another copy (e.g. the notebook's cache/prep_ALL.qs2).
# Output: FIG2_OUT/panel_B_umap_by_type_bare.{pdf,png,svg} (figure panel; also _songlabels with titles), panel_B_umap_coordinates.csv.gz,
#         FIG2_OUT/panel_B_reciprocal_score_bare.{pdf,png,svg} (figure panel; also panel_B_reciprocal_score with its colour key),
#         panel_B_reciprocal_scores.csv
#
#   Rscript figures/figure2/panel_B_song_nonsong_integration.R
source(here::here("figures/figure2/common.R"))
suppressMessages({library(Seurat); library(qs2)})
options(future.globals.maxSize = Inf)
set.seed(10)

LABEL_COL <- "celltype_hybrid"
SONG_PATTERN <- "HVC|LMAN|-RA$"   # song-nucleus types; every other glutamatergic type is non-song
SET_PATTERN <- "^Glut-"           # the "ALL" analysis: every glutamatergic type
N_FEATURES <- 3000
DIMS <- 1:30
REDUCTION <- "cca"
K_ANCHOR <- 20
K_WEIGHT <- 50                    # TransferData
PREP_FILE <- Sys.getenv("FIG2_B_PREP", file.path(FIG2_OUT, "panel_B_prep_ALL.qs2"))

n_threads <- min(16L, max(1L, parallel::detectCores() - 4L))
RhpcBLASctl::blas_set_num_threads(n_threads)
RcppParallel::setThreadOptions(numThreads = n_threads)

# ---- per-dataset SCTransform (song_nonsong_integration_hybrid.qmd) ----------------------------------------------------------------------------------
if (!file.exists(PREP_FILE)) {
  obj <- qs_read(GLUT_OBJ, nthreads = 8)
  obj$label <- droplevels(as.factor(obj@meta.data[[LABEL_COL]]))
  obj$set <- ifelse(grepl(SONG_PATTERN, obj$label), "song", "nonsong")
  obj <- subset(obj, cells = colnames(obj)[grepl(SET_PATTERN, obj$label)])
  obj$label <- droplevels(obj$label)
  # the RNA assay is BPCells-backed on disk; materialize the counts into a plain object
  cnt <- as(LayerData(obj, assay = "RNA", layer = "counts"), "dgCMatrix")
  mem <- CreateSeuratObject(cnt, meta.data = obj@meta.data[, c("set", "label", "sample_name", "position"), drop = FALSE])
  lst <- SplitObject(mem, split.by = "set")[c("song", "nonsong")]
  lst <- lapply(lst, function(x) SCTransform(x, verbose = FALSE))
  features <- SelectIntegrationFeatures(lst, nfeatures = N_FEATURES, verbose = FALSE)
  lst <- PrepSCTIntegration(lst, anchor.features = features, verbose = FALSE)
  lst <- lapply(lst, function(x) RunPCA(x, features = features, npcs = max(DIMS), verbose = FALSE))
  prep <- list(list = lst, features = features,
               labels = setNames(as.character(obj$label), colnames(obj)), libs = setNames(as.character(obj$sample_name), colnames(obj)))
  qs_save(prep, PREP_FILE, nthreads = 8)
  rm(obj, cnt, mem, lst); invisible(gc())
}
prep <- qs_read(PREP_FILE, nthreads = 8)
lab <- prep$labels
message("cells: ", length(lab), " (", sum(grepl(SONG_PATTERN, lab)), " song), types: ", length(unique(lab)))

# ---- UMAPs (song_nonsong_umap_hybrid.qmd) -----------------------------------------------------------------------------------------------------------
umap_tbl <- function(o, embedding) {
  as.data.frame(Embeddings(o, "umap")) %>% rownames_to_column("cell") %>% setNames(c("cell", "umap_1", "umap_2")) %>% mutate(embedding = embedding)
}
mg <- merge(prep$list[[1]], prep$list[[2]])
DefaultAssay(mg) <- "SCT"
VariableFeatures(mg) <- prep$features
mg <- mg %>% ScaleData(features = prep$features, verbose = FALSE) %>%
  RunPCA(features = prep$features, npcs = max(DIMS), verbose = FALSE) %>% RunUMAP(dims = DIMS, verbose = FALSE)
umap_raw <- umap_tbl(mg, "uncorrected")
rm(mg); invisible(gc())

anchors <- FindIntegrationAnchors(prep$list, anchor.features = prep$features, normalization.method = "SCT",
                                  reduction = REDUCTION, dims = DIMS, k.anchor = K_ANCHOR, verbose = FALSE)
k_weight_int <- min(100, min(vapply(prep$list, ncol, integer(1))) - 1)
int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = DIMS, k.weight = k_weight_int, verbose = FALSE)
DefaultAssay(int) <- "integrated"
int <- int %>% RunPCA(npcs = max(DIMS), verbose = FALSE) %>% RunUMAP(dims = DIMS, verbose = FALSE)
umap_int <- umap_tbl(int, "integrated")
rm(int, anchors); invisible(gc())

umap <- bind_rows(umap_raw, umap_int) %>%
  mutate(label = unname(lab[cell]), dataset = if_else(grepl(SONG_PATTERN, label), "song", "non-song"), library = unname(prep$libs[cell]),
         embedding = factor(embedding, levels = c("uncorrected", "integrated")))
write_csv(umap, file.path(FIG2_OUT, "panel_B_umap_coordinates.csv.gz"))

# ---- reciprocal label transfer (song_nonsong_reciprocal_transfer_hybrid.qmd) ------------------------------------------------------------------------
#' frac x mean_score matrix of one transfer direction: query types as rows, reference types as columns
transfer_matrix <- function(ref, query, query_labels) {
  anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", features = prep$features,
                                 reduction = REDUCTION, dims = DIMS, k.anchor = K_ANCHOR, verbose = FALSE)
  n_anchor_cells <- min(n_distinct(anchors@anchors[, "cell1"]), n_distinct(anchors@anchors[, "cell2"]))
  pred <- TransferData(anchorset = anchors, refdata = ref$celltype, weight.reduction = REDUCTION, dims = DIMS,
                       k.weight = min(K_WEIGHT, n_anchor_cells - 1), verbose = FALSE)
  score_cols <- setdiff(grep("^prediction\\.score\\.", colnames(pred), value = TRUE), "prediction.score.max")
  scores <- pred[, score_cols, drop = FALSE]
  ref_types <- sort(unique(ref$celltype))
  colnames(scores) <- ref_types[match(gsub("^prediction\\.score\\.", "", score_cols), make.names(ref_types))]
  stopifnot(!any(is.na(colnames(scores))))
  as_tibble(scores, rownames = "cell") %>%
    mutate(query_type = unname(query_labels[cell]), predicted = pred$predicted.id) %>%
    pivot_longer(-c(cell, query_type, predicted), names_to = "ref_type", values_to = "score") %>%
    group_by(query_type, ref_type) %>%
    summarise(frac = mean(predicted == ref_type), mean_score = mean(score), .groups = "drop") %>%
    mutate(value = frac * replace_na(mean_score, 0)) %>%
    select(query_type, ref_type, value) %>%
    pivot_wider(names_from = ref_type, values_from = value) %>% column_to_rownames("query_type") %>% as.matrix()
}
song <- prep$list[["song"]]; song$celltype <- unname(lab[colnames(song)])
nonsong <- prep$list[["nonsong"]]; nonsong$celltype <- unname(lab[colnames(nonsong)])
fwd <- transfer_matrix(ref = nonsong, query = song, setNames(song$celltype, colnames(song)))         # song cells labelled from non-song
rev <- transfer_matrix(ref = song, query = nonsong, setNames(nonsong$celltype, colnames(nonsong)))  # non-song cells labelled from song
song_types <- rownames(fwd); nonsong_types <- colnames(fwd)
Fm <- fwd[song_types, nonsong_types, drop = FALSE]
Rm <- matrix(0, length(song_types), length(nonsong_types), dimnames = list(song_types, nonsong_types))
shared_r <- intersect(nonsong_types, rownames(rev)); shared_c <- intersect(song_types, colnames(rev))
Rm[shared_c, shared_r] <- t(rev[shared_r, shared_c, drop = FALSE])
recip <- tibble(song_type = rep(song_types, times = length(nonsong_types)), nonsong_type = rep(nonsong_types, each = length(song_types)),
                fwd = as.vector(Fm), rev = as.vector(Rm), reciprocal = sqrt(as.vector(Fm) * as.vector(Rm)))
write_csv(recip, file.path(FIG2_OUT, "panel_B_reciprocal_scores.csv"))
rm(song, nonsong, prep); invisible(gc())

# ---- figures ----------------------------------------------------------------------------------------------------------------------------------------
# The fixed-size panels below are laid out with text measured on the open device. The notebooks drew them on knitr's ragg device (200 dpi);
# open the same here, otherwise R's default pdf() device measures Arial with substitute metrics and the heatmap comes out smaller.
ragg::agg_png(tempfile(fileext = ".png"), width = 2.5, height = 2.5, units = "in", res = 200)
# Song types drop the division prefix; non-song types are abbreviated (D2-n, C1H-n)
short <- function(x) if_else(grepl(SONG_PATTERN, x), str_replace(x, "^Glut-(DACH2|CACNA1H)-", ""),
                             x %>% str_replace("^Glut-DACH2-", "D2-") %>% str_replace("^Glut-CACNA1H-", "C1H-") %>% str_replace("^Glut-", ""))
song_position <- function(x) case_when(grepl("LMAN", x) ~ "lman", grepl("HVC", x) ~ "hvc", grepl("RA", x) ~ "ra", TRUE ~ NA_character_)

#' PDF (cairo), PNG and SVG (one <text> element per label) of a gtable, at its own fixed size
fig_save_gt <- function(gt, path_noext, dpi = 600) {
  fig_check_font()
  w <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
  h <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
  ggsave(paste0(path_noext, ".pdf"), gt, width = w, height = h, device = cairo_pdf, family = FIG_FONT)
  ggsave(paste0(path_noext, ".png"), gt, width = w, height = h, dpi = dpi, bg = FIG_SURFACE, device = ragg::agg_png)
  ggsave(paste0(path_noext, ".svg"), gt, width = w, height = h, device = svglite::svglite, bg = FIG_SURFACE)
  message("saved ", path_noext)
}

# UMAPs: 1.4 in square panels, points rasterized at 600 dpi, song types in their nucleus colour, non-song grey
umap_plot <- function(d, bare) {
  u <- unique(d$label)
  cols <- setNames(if_else(grepl(SONG_PATTERN, u), unname(POSITION_COLORS[song_position(u)]), FIG_INK_MUTED), u)
  lab_d <- d %>% group_by(embedding) %>%
    group_modify(~ .x %>% group_by(label) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop") %>% mutate(lab = short(label))) %>%
    ungroup() %>% filter(grepl(SONG_PATTERN, label))
  p <- ggplot(d, aes(umap_1, umap_2)) +
    ggrastr::rasterise(geom_point(aes(colour = label), size = 0.12, stroke = 0, alpha = 0.55, show.legend = FALSE), dpi = 600, dev = "ragg") +
    ggrepel::geom_text_repel(data = lab_d, aes(label = lab), size = fig_pt(5), family = FIG_FONT, colour = FIG_INK_PRIMARY, segment.size = 0.15,
                             segment.colour = FIG_INK_MUTED, min.segment.length = 0.1, max.overlaps = Inf, box.padding = 0.1, seed = 1) +
    facet_wrap(~embedding, nrow = 1) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
    theme_fig() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY))
  p <- if (bare) p + labs(x = NULL, y = NULL) + theme(strip.text = element_blank(), plot.margin = margin(1, 1, 1, 1))
       else p + labs(x = "UMAP 1", y = "UMAP 2", title = "ALL glutamatergic types",
                     subtitle = sprintf("%s, k.anchor = %d | song types coloured, non-song grey", toupper(REDUCTION), K_ANCHOR))
  egg::set_panel_size(p, width = unit(1.4, "in"), height = unit(1.4, "in"))
}
fig_save_gt(umap_plot(umap, bare = TRUE), file.path(FIG2_OUT, "panel_B_umap_by_type_bare"))
fig_save_gt(umap_plot(umap, bare = FALSE), file.path(FIG2_OUT, "panel_B_umap_by_type_songlabels"))

# Heatmap: 0.1 in cells, Blues ramp from 0 to the maximum; song columns HVC, LMAN, RA; non-song rows D2-1..D2-8, C1H-1..C1H-4
blues <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"), space = "rgb")(256)
score_panel <- function(d, bare) {
  u <- unique(d$song_type); lv_song <- rev(u[order(match(song_position(u), c("hvc", "lman", "ra")), u)])
  u <- unique(d$nonsong_type)
  lv_ns <- u[order(match(str_extract(u, "CACNA1H|DACH2"), c("DACH2", "CACNA1H")), as.integer(str_extract(u, "[0-9]+$")))]
  d <- d %>% mutate(song_type = factor(song_type, levels = lv_song), nonsong_type = factor(nonsong_type, levels = lv_ns))
  p <- ggplot(d, aes(song_type, nonsong_type, fill = reciprocal)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_gradientn(colours = blues, limits = c(0, max(d$reciprocal)), name = "score",
                         guide = guide_colourbar(barwidth = unit(0.06, "in"), barheight = unit(0.45, "in"), title.position = "top")) +
    scale_x_discrete(labels = short, limits = rev(lv_song), expand = fig_expand(0)) +
    scale_y_discrete(labels = short, limits = rev(lv_ns), expand = fig_expand(0)) +
    labs(x = NULL, y = NULL) +
    theme_fig() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = rev(unname(POSITION_COLORS[song_position(lv_song)]))),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = if (bare) "none" else "right",
          legend.margin = margin(0, 0, 0, 0), legend.title = element_text(size = FIG_PT_AXIS_TEXT), plot.margin = margin(1, 1, 1, 1))
  egg::set_panel_size(p, width = unit(length(lv_song) * 0.1, "in"), height = unit(length(lv_ns) * 0.1, "in"))
}
fig_save_gt(score_panel(recip, bare = TRUE), file.path(FIG2_OUT, "panel_B_reciprocal_score_bare"))
fig_save_gt(score_panel(recip, bare = FALSE), file.path(FIG2_OUT, "panel_B_reciprocal_score"))
invisible(dev.off())
