#!/usr/bin/env Rscript
# Figure 2, panel A -- hierarchical clustering of the nidopallial (DACH2) and arcopallial (CACNA1H) glutamatergic cell types after removing song DEGs.
#
# Glutamatergic cells with the hybrid labels, without the developmental types (Glut-NSC, Glut-NB, Glut-Im), Glut-DACH2-HVCra-Int, Glut-SATB2-1 and
# Glut-GABA. Within each division (DACH2 = nidopallial dissections, CACNA1H = arcopallial), markers of every cell type are called on pseudobulk
# samples (cell type x bird; Seurat FindAllMarkers, DESeq2, SCT counts). Genes called for a song type (HVCra, HVCx, LMANco, LMANsh; RA) at adjusted
# p < 0.1 are removed; of the rest, the 100 strongest up- and 100 strongest down-regulated markers (adjusted p < 0.1) of each type are the features.
# Each feature's specificity across the division's types (calc_gene_specificity on the log1p of summed SCT counts per type) is clustered by Ward's
# method (ward.D) on correlation distance; node labels are pvclust approximately-unbiased support (%, 1,000 multiscale bootstrap replicates). Tips
# are coloured by the dissected region of the cell type's original (pre-hybrid) cluster name.
# Same computation and drawing as snrna/trees/celltypes_hclust_glut_robustness_hybrid.qmd at its manuscript threshold (0.1), which wrote it as page
# 4 of correlation_ward.D/trees_nido_thresholds.pdf and trees_arco_thresholds.pdf.
#
# Bootstrap: the notebook ran pvclust in parallel without a seed, so its support values cannot be regenerated exactly; this script seeds the
# bootstrap (iseed = 1). To draw the support values of the run in the figure instead, set FIG2_PVCLUST_SWEEP to that run's saved pvclust objects
# (correlation_ward.D/pvclust_sweep.qs2); the script then checks that the saved tree equals the recomputed one before using them.
#
# Input: the snRNA-seq object with the hybrid labels (FIG2_SNRNA_OBJ; snrna/naming/hybrid_division_naming.qmd).
# Outputs: FIG2_OUT/panel_A_tree_nido.pdf, panel_A_tree_arco.pdf; panel_A_tree_support.csv (per-node support), panel_A_features.csv (genes clustered)
#
#   Rscript figures/figure2/panel_A_glut_trees.R
#   FIG2_PVCLUST_SWEEP=snrna/trees/celltypes_hclust_glut_robustness_hybrid/correlation_ward.D/pvclust_sweep.qs2 Rscript figures/figure2/panel_A_glut_trees.R
source(here::here("figures/figure2/common.R"))
suppressMessages({library(Seurat); library(qs2); library(pvclust); library(ggtree); library(ape)})

SNRNA_OBJ <- Sys.getenv("FIG2_SNRNA_OBJ", here::here("snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2"))
PVCLUST_SWEEP <- Sys.getenv("FIG2_PVCLUST_SWEEP", "")

ASSAY <- "SCT"
NSIG <- 100               # markers kept per cell type and direction
FEATURE_THRESH <- 0.1     # adjusted p for a marker to be a clustering feature
DEG_THRESH <- 0.1         # adjusted p for a song-type marker to be removed (the manuscript threshold)
NBOOT <- 1000
N_CORES <- min(20L, max(1L, parallel::detectCores() - 4L))
GLUT_EXCLUDE <- c("Glut-NSC", "Glut-NB", "Glut-Im", "Glut-DACH2-HVCra-Int", "Glut-SATB2-1", "Glut-GABA")
REGIONS <- list(NIDO = list(cols = "^Glut-DACH2-", song = "HVC|LMAN"),
                ARCO = list(cols = "^Glut-CACNA1H-", song = "^Glut-CACNA1H-RA$"))

#' Gene specificity across cell types (/opt/colquittlab/utils/R/stats.R)
calc_gene_specificity <- function(x, scale = FALSE) {
  if (scale) x <- (x - min(x)) / (max(x) - min(x))
  x_norm <- x / sum(x)
  x_norm * log(x_norm / mean(x_norm) + .01)
}

# ---- cells ------------------------------------------------------------------------------------------------------------------------------------------
obj <- qs_read(SNRNA_OBJ, nthreads = 8)
ct <- obj$celltype_hybrid
obj <- subset(obj, cells = colnames(obj)[!is.na(ct) & grepl("^Glut", ct) & !(ct %in% GLUT_EXCLUDE)])
stopifnot("doublet assignments present" = !any(grepl("/", obj$assignment)))
obj$region <- case_when(obj$position %in% c("ra", "arco") ~ "ARCO", obj$position %in% c("lman", "nr", "hvc", "nc") ~ "NIDO")
obj$orig_cluster <- as.character(obj$cluster)
obj$cluster <- droplevels(as.factor(obj$celltype_hybrid))
# pseudobulk replicate = bird; souporcell's assignment only identifies a bird within its library
obj$replicate_lib <- paste(obj$sample_name, obj$assignment, sep = "_")

# division prefix and dissection agree for every type (DACH2 -> NIDO, CACNA1H -> ARCO)
maj <- obj@meta.data %>% count(cluster, region) %>% group_by(cluster) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()
stopifnot(all(maj$region == if_else(grepl("^Glut-DACH2-", maj$cluster), "NIDO", "ARCO")))

# tip colour: dissected region parsed from the original cluster name (Glut-<Region>-n; Nido counts as NC)
tip_position <- obj@meta.data %>% distinct(cluster, orig_cluster) %>%
  mutate(position = str_extract(orig_cluster, "Glut-([A-Za-z]+)", group = 1),
         position = tolower(if_else(position == "Nido", "NC", position))) %>%
  distinct(cluster, position)
stopifnot(!anyDuplicated(tip_position$cluster), all(tip_position$position %in% names(POSITION_COLORS)))

# ---- markers and per-type expression ----------------------------------------------------------------------------------------------------------------
markers <- map(names(REGIONS), function(reg) {
  agg <- AggregateExpression(subset(obj, region == reg), assay = ASSAY, group.by = c("cluster", "replicate_lib"), return.seurat = TRUE)
  Idents(agg) <- agg$cluster
  FindAllMarkers(agg, assay = ASSAY, test.use = "DESeq2", only.pos = FALSE, verbose = FALSE)
}) %>% set_names(names(REGIONS))
Idents(obj) <- obj$cluster
avg <- AggregateExpression(obj, assays = ASSAY)[[ASSAY]]
avg <- log1p(avg)
rm(obj); invisible(gc())

# ---- features, specificity, tree --------------------------------------------------------------------------------------------------------------------
features <- map(names(REGIONS), function(reg) {
  m <- markers[[reg]]
  removed <- m %>% filter(grepl(REGIONS[[reg]]$song, cluster), p_val_adj < DEG_THRESH) %>% pull(gene) %>% unique()
  genes <- m %>% filter(!(gene %in% removed)) %>% mutate(sign = avg_log2FC > 0) %>% filter(p_val_adj < FEATURE_THRESH) %>%
    group_by(cluster, sign) %>% top_n(NSIG, abs(avg_log2FC)) %>% ungroup() %>% distinct(gene, .keep_all = TRUE) %>% pull(gene)
  list(removed = removed, genes = genes)
}) %>% set_names(names(REGIONS))
imap_dfr(features, ~ tibble(region = .y, gene = .x$genes)) %>% write_csv(file.path(FIG2_OUT, "panel_A_features.csv"))

spec_matrix <- function(genes, col_pattern) {
  m <- as.matrix(avg[genes, grepl(col_pattern, colnames(avg)), drop = FALSE])
  m <- t(apply(m, 1, calc_gene_specificity, scale = TRUE))
  m[is.nan(m)] <- 0
  m
}
saved <- if (nzchar(PVCLUST_SWEEP)) qs_read(PVCLUST_SWEEP) else NULL
pv <- map(names(REGIONS), function(reg) {
  mat <- spec_matrix(features[[reg]]$genes, REGIONS[[reg]]$cols)
  p <- pvclust(mat, method.dist = "correlation", method.hclust = "ward.D", nboot = NBOOT, parallel = as.integer(N_CORES), iseed = 1, quiet = TRUE)
  if (!is.null(saved)) {
    s <- saved[[paste(reg, format(DEG_THRESH, scientific = FALSE), sep = "_")]]$pv
    stopifnot("saved tree differs from the recomputed one" =
                identical(s$hclust$merge, p$hclust$merge), identical(s$hclust$labels, p$hclust$labels),
              isTRUE(all.equal(s$hclust$height, p$hclust$height)))
    p <- s
  }
  p
}) %>% set_names(names(REGIONS))

# ---- drawing (as in the notebook's support_tree()) --------------------------------------------------------------------------------------------------
hclust_clades <- function(hc) {   # tip sets of the internal nodes, in hclust merge order (= pvclust's $edges rows)
  out <- vector("list", nrow(hc$merge))
  for (i in seq_len(nrow(hc$merge))) {
    tips <- character(0)
    for (j in hc$merge[i, ]) tips <- c(tips, if (j < 0) hc$labels[-j] else out[[j]])
    out[[i]] <- sort(tips)
  }
  out
}
clade_key <- function(tips) paste(sort(tips), collapse = "|")
phylo_tipsets <- function(phy) {
  desc <- vector("list", max(phy$edge))
  for (i in seq_along(phy$tip.label)) desc[[i]] <- phy$tip.label[i]
  for (e in rev(seq_len(nrow(phy$edge)))) desc[[phy$edge[e, 1]]] <- c(desc[[phy$edge[e, 1]]], desc[[phy$edge[e, 2]]])
  desc
}
edge_values_to_nodes <- function(hc, values) {
  ts <- phylo_tipsets(ape::as.phylo(hc))
  tibble(node = seq_along(ts), key = map_chr(ts, ~ if (is.null(.x)) NA_character_ else clade_key(.x))) %>%
    filter(!is.na(key)) %>%
    left_join(tibble(key = map_chr(hclust_clades(hc), clade_key), value = values), by = "key") %>%
    filter(!is.na(value))
}

TIP_PT <- 6; SUPPORT_PT <- 5
TIP_SCALE <- TIP_PT * 1.45 / 72   # panel inches per tip
PANEL_W <- 2
support_tree <- function(hc, values, title, subtitle) {
  tips <- data.frame(labels = hc$labels) %>% left_join(tip_position, by = c("labels" = "cluster")) %>%
    mutate(col = unname(POSITION_COLORS[position]))
  nodes <- edge_values_to_nodes(hc, values) %>% transmute(node, support_lab = as.character(round(value * 100)))
  g <- ggtree(hc, branch.length = "none", linewidth = 0.5 / .pt) %<+% tips %<+% nodes +
    geom_tiplab(aes(color = col), angle = 0, hjust = 0, offset = 0, show.legend = FALSE, size = TIP_PT / .pt) +
    theme_dendrogram(plot.margin = margin(20, 20, 20, 20)) +
    scale_color_identity() +
    hexpand(5) +
    labs(title = title, subtitle = subtitle) +
    theme(plot.title = element_text(size = 9, face = "plain"), plot.subtitle = element_text(size = 7)) +
    geom_text2(aes(label = support_lab, subset = !isTip), hjust = 1.15, vjust = -0.4, size = SUPPORT_PT / .pt, color = "grey30") +
    ggexpand(0.25, direction = -1)
  egg::set_panel_size(g, width = unit(PANEL_W, "in"), height = unit(length(hc$labels) * TIP_SCALE, "in"))
}

support <- list()
for (reg in names(REGIONS)) {
  p <- pv[[reg]]
  gt <- support_tree(p$hclust, p$edges$au, title = sprintf("%s\nsong DEGs removed at p_adj < %s", reg, format(DEG_THRESH, scientific = FALSE)),
                     subtitle = sprintf("%d genes removed, %d genes clustered\nbranch labels = pvclust AU (%%)",
                                        length(features[[reg]]$removed), length(features[[reg]]$genes)))
  out <- file.path(FIG2_OUT, sprintf("panel_A_tree_%s.pdf", tolower(reg)))
  pdf(out, height = length(p$hclust$labels) * TIP_SCALE + 1.4, width = PANEL_W + 2.0)
  grid::grid.draw(gt)
  invisible(dev.off())
  support[[reg]] <- tibble(region = reg, clade = map_chr(hclust_clades(p$hclust), clade_key), au = p$edges$au, bp = p$edges$bp)
  message("saved ", out)
}
bind_rows(support) %>% write_csv(file.path(FIG2_OUT, "panel_A_tree_support.csv"))
