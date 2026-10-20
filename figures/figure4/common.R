# figures/figure4/common.R -- inputs and shared helpers for the Figure 4 panel scripts.
#
# Every panel script sources this file (and the repository's config/ files through it) and nothing else in this folder, so each panel can be
# regenerated on its own. Inputs are the outputs of the upstream pipelines named in README.md; edit the paths below (or set the environment variables)
# to point at them. Panels are written to FIG4_OUT (default: figures/figure4/output/, not tracked).
suppressMessages({library(here); library(tidyverse)})
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))

FIG4_OUT <- Sys.getenv("FIG4_OUT", here::here("figures/figure4/output"))
dir.create(FIG4_OUT, showWarnings = FALSE, recursive = TRUE)

# ---- inputs -----------------------------------------------------------------------------------------------------------------------------------------
# SCENIC+ run on the hybrid-labeled multiome data (all cell types), config37 (grn/ra-arco-hvc-nc/scenicplus/)
HYBRID_RUN <- file.path(SCENICPLUS_JUPYTER_DIR, "ra-arco-hvc-nc_hybrid")
SCPLUS_CONFIG <- "config37"
SCPLUS_OUTS <- file.path(HYBRID_RUN, "scenicplus", SCPLUS_CONFIG, "outs")
SCPLUS_H5MU <- Sys.getenv("FIG4_SCPLUS_H5MU", if (file.exists(file.path(SCPLUS_OUTS, "scplusmdata.h5mu"))) file.path(SCPLUS_OUTS, "scplusmdata.h5mu") else file.path(SCPLUS_OUTS, "scplusmdata_slim.h5mu"))
# ArchR project copy on the SCENIC+ consensus regions, with the saved marker tests (grn/ra-arco-hvc-nc/pycisTopic/make_song_pair_dars_archr.R and
# make_interneuron_dars_archr.R write markers_<FG>_VS_<BG>.rds here)
ARCHR_CONSENSUS <- Sys.getenv("FIG4_ARCHR_CONSENSUS", file.path(HYBRID_RUN, "archr_consensus"))
# Hybrid-labeled multiome Seurat object with the UMAP reductions (multiome/seurat/reduction_viz/combined_all_umap_hybrid.R)
SEURAT_HYBRID_OBJ <- Sys.getenv("FIG4_SEURAT_OBJ", here::here("multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"))
UMAP_REDUCTION <- "dims25nn30mindist0.3"

# ---- cell types -------------------------------------------------------------------------------------------------------------------------------------
RA <- "Glut-CACNA1H-RA"; C1H1 <- "Glut-CACNA1H-1"; PV1 <- "GABA-MGE-PVALB-1"; LGE1 <- "GABA-LGE-1"
# column order of the cell-type heatmaps (as in grn/analysis/scenicplus_hybrid_all_template.qmd)
celltype_order <- function(cl) unique(c(sort(grep("HVC", cl, value = TRUE)), sort(grep("DACH2-[0-9]", cl, value = TRUE)), sort(grep("CACNA1H-RA", cl, value = TRUE)),
                                        sort(grep("CACNA1H-[0-9]", cl, value = TRUE)), intersect(c("Glut-Im", "Glut-NB", "Glut-NSC"), cl), sort(setdiff(grep("GABA", cl, value = TRUE), "GABA-Im")),
                                        intersect("GABA-Im", cl), sort(grep("Astro", cl, value = TRUE)), sort(grep("Epen", cl, value = TRUE)), sort(grep("Oligo", cl, value = TRUE)),
                                        sort(grep("OPC", cl, value = TRUE)), sort(grep("Micro", cl, value = TRUE)), sort(grep("Endo", cl, value = TRUE))))
# LOC genes named by protein homology, or ncRNA-N (grn/analysis/loc_gene_labels_lonStrDom2.csv)
LABEL_MAP <- c(LOC110468166 = "NTM", LOC110470685 = "TLE4", LOC110474504 = "VIPR1", LOC116183613 = "ncRNA-1", LOC110482360 = "LAMA3")
gene_label <- function(g) ifelse(g %in% names(LABEL_MAP), LABEL_MAP[g], g)

# ---- SCENIC+ outputs (read straight from the h5mu) -------------------------------------------------------------------------------------------------
h5_factor <- function(x) factor(x$codes, levels = 0:(length(x$categories) - 1), labels = x$categories)

#' eRegulon table (direct + extended), as written by SCENIC+ into the h5mu; `signs` keeps the +/+ (and/or -/-) eRegulons
read_eregulons <- function(h5mu = SCPLUS_H5MU, signs = "+/+") {
  suppressMessages(library(rhdf5))
  e <- map_dfr(c("direct", "extended"), function(t) as_tibble(do.call(cbind.data.frame, h5read(h5mu, paste0("/uns/", t, "_e_regulon_metadata")))))
  e %>% filter(sapply(eRegulon_name, function(n) any(endsWith(n, signs))))
}

#' Cell x eRegulon gene-based AUC (direct and extended), with each cell's cluster, library (sample_short) and bird (assignment)
read_gene_auc <- function(h5mu = SCPLUS_H5MU) {
  suppressMessages(library(rhdf5))
  one <- function(t) { cur <- h5read(h5mu, paste0("/mod/", t, "_gene_based_AUC")); m <- data.frame(t(cur[["X"]]), check.names = FALSE); colnames(m) <- cur[["var"]]$`_index`; rownames(m) <- cur[["obs"]]$Cell; m }
  auc <- bind_cols(one("direct"), one("extended")) %>% rownames_to_column("_index")
  obs <- h5read(h5mu, "/mod/scRNA_counts/obs")
  grp <- tibble(`_index` = obs$`_index`, cluster = h5_factor(obs$cluster), sample_short = h5_factor(obs$sample_short), assignment = h5_factor(obs$assignment))
  auc %>% left_join(grp, by = "_index") %>% filter(!is.na(cluster))
}

#' Per-TF positive (+/+) regulon AUC in long form: the direct and extended regulons of a TF share the name (their AUC values are pooled), as in panel A
positive_auc_long <- function(auc, by = "cluster") {
  pos <- auc %>% select(all_of(by), matches("\\+\\/\\+"))
  colnames(pos) <- sub("_(direct|extended)_\\+\\/\\+_\\([0-9]+[a-z]\\)", "", colnames(pos))
  pos %>% pivot_longer(-all_of(by), names_to = "tf")
}

# ---- ArchR differential accessibility ---------------------------------------------------------------------------------------------------------------
#' Region log2 fold changes of a saved ArchR marker test (getMarkerFeatures, bias-matched Wilcoxon); regions named chr:start-end in BED coordinates
archr_lfc <- function(fg, bg) {
  suppressMessages(library(SummarizedExperiment))
  se <- readRDS(file.path(ARCHR_CONSENSUS, paste0("markers_", fg, "_VS_", bg, ".rds")))
  rr <- rowData(se)
  setNames(as.numeric(assay(se, "Log2FC")[, 1]), paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end))
}

# ---- pseudobulk RNA differential expression ---------------------------------------------------------------------------------------------------------
#' Pseudobulk counts: one sample per library x bird x cell type (sample_short x assignment x cluster)
pseudobulk_counts <- function(obj) {
  suppressMessages({library(Seurat); library(Matrix)})
  meta <- obj@meta.data; meta$cluster <- meta$cluster_hybrid
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
  f <- factor(paste(meta$sample_short, meta$assignment, meta$cluster, sep = "|"))
  pb <- cnt %*% sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  colnames(pb) <- levels(f)
  parts <- do.call(rbind, strsplit(levels(f), "|", fixed = TRUE))
  list(pb = pb, info = tibble(id = levels(f), lib = parts[, 1], cluster = parts[, 3], n = as.integer(table(f))))
}

#' DESeq2 on pseudobulk samples with >= min_cells cells: fg cluster(s) vs bg cluster(s), design ~ library + group when estimable (else ~ group)
pseudobulk_de <- function(pbc, fg, bg, min_cells = 10) {
  suppressMessages(library(DESeq2))
  s <- pbc$info %>% filter(cluster %in% c(fg, bg), n >= min_cells) %>% mutate(g = ifelse(cluster %in% fg, "fg", "bg"))
  cd <- data.frame(row.names = s$id, lib = factor(s$lib), grp = factor(s$g, levels = c("bg", "fg")))
  stopifnot(sum(cd$grp == "fg") >= 2, sum(cd$grp == "bg") >= 2)
  mm <- model.matrix(~ lib + grp, cd)
  design <- if (qr(mm)$rank == ncol(mm)) ~ lib + grp else ~ grp
  dds <- DESeqDataSetFromMatrix(round(as.matrix(pbc$pb[, s$id])), cd, design)
  dds <- DESeq(dds[rowSums(counts(dds) >= 5) >= 2, ], quiet = TRUE)
  r <- results(dds, contrast = c("grp", "fg", "bg"))
  message(sprintf("  %s vs %s: %d vs %d pseudobulk samples, design %s, %d genes", paste(fg, collapse = "+"), paste(bg, collapse = "+"), sum(cd$grp == "fg"), sum(cd$grp == "bg"), format(design), nrow(r)))
  tibble(gene = rownames(r), lfc = r$log2FoldChange, padj = r$padj)
}

#' The two MAFB on/off RNA contrasts of panels F and G for every gene: RA vs C1H-1 and PVALB-1 vs LGE-1
mafb_onoff_rna_contrasts <- function() {
  suppressMessages(library(qs2))
  obj <- qs_read(SEURAT_HYBRID_OBJ, nthreads = 8)
  pbc <- pseudobulk_counts(obj); rm(obj); invisible(gc())
  pseudobulk_de(pbc, RA, C1H1) %>% select(gene, lfc_RA = lfc, padj_RA = padj) %>%
    full_join(pseudobulk_de(pbc, PV1, LGE1) %>% select(gene, lfc_PV1 = lfc, padj_PV1 = padj), by = "gene")
}

# ---- output ------------------------------------------------------------------------------------------------------------------------------------------
#' Write a ggplot (optionally with a fixed panel size) as PDF, transparent PNG and SVG (one <text> element per label)
save_gg <- function(p, name, width, height) {
  fig_save(p, file.path(FIG4_OUT, name), width = width, height = height)
  ggsave(file.path(FIG4_OUT, paste0(name, ".png")), p, width = width, height = height, dpi = 600, bg = "transparent")
  ggsave(file.path(FIG4_OUT, paste0(name, ".svg")), p, width = width, height = height, bg = "transparent",
         device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
  message("saved ", file.path(FIG4_OUT, name))
}
