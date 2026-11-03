# overexpression/reanalysis/oe_helpers.R -- shared helpers for the MAFB/EMX2 overexpression reanalysis.
#
#   source(here::here("overexpression/reanalysis/oe_helpers.R"))
#
# Defines values and functions only. Two groups of helpers:
#   - pseudobulk DESeq2, following figures/figure4/common.R (pseudobulk_counts / pseudobulk_de): one replicate per
#     library x individual x group, >= MIN_CELLS cells, rows kept when >= 5 counts in >= 2 samples;
#   - gene-signature scoring with an expression-matched null, copied from the songbird_dev repo
#     (/ssd/shared/multiome/dev_full1/lonStrDom3/seurat/dev-full1_cellbender.0.1_multiome_seurat_preprocess/nido/,
#     calibrated_annotation_helpers.R and core_signature_helpers.R, merged to main in PR #5). They are copied rather than
#     sourced so this analysis does not break when that repo moves; the function bodies are unchanged except where noted.
suppressMessages({library(tidyverse); library(Matrix)})
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))

# UCell and AUCell are not in this project's renv library. The dev_full1 project's R 4.6 library has both (the versions
# the scoring method was developed with); append it, so this project's own packages still take precedence.
DEV_FULL1_RLIB <- Sys.getenv("DEV_FULL1_RLIB",
                             "/ssd/shared/multiome/dev_full1/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
if (dir.exists(DEV_FULL1_RLIB)) .libPaths(unique(c(.libPaths(), DEV_FULL1_RLIB)))

OE_DIR <- here::here("overexpression/reanalysis")
# Large objects and figures go here (gitignored by extension); small tables go to OE_DIR/tables (tracked).
OE_OUT <- Sys.getenv("OE_OUT", file.path(OE_DIR, "output"))
OE_TABLES <- file.path(OE_DIR, "tables")
for (d in c(OE_OUT, OE_TABLES)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Adult reference (zebra/Bengalese finch snRNA, hybrid labels) -- written by snrna/naming/hybrid_division_naming.qmd.
# Gitignored, so in a worktree fall back to the main checkout's copy.
ref_obj_path <- function() {
  rel <- "snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2"
  cands <- here::here(rel)
  main <- tryCatch(system2("git", c("-C", shQuote(here::here()), "rev-parse", "--path-format=absolute", "--git-common-dir"),
                           stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (length(main) == 1 && nzchar(main)) cands <- c(cands, file.path(dirname(main), rel))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else cands[1]
}

# Song-neuron axes: song cluster vs its non-song sister.
AXES <- tribble(
  ~axis,   ~song,               ~sister,
  "RA",    "Glut-CACNA1H-RA",   "Glut-CACNA1H-1",
  "HVCra", "Glut-DACH2-HVCra",  "Glut-DACH2-1",
  "HVCx",  "Glut-DACH2-HVCx",   "Glut-DACH2-4"
)
# Adult glutamatergic clusters for the specificity gate: every Glut-* label except the immature/progenitor ones and the
# mixed Glut-GABA cluster.
ADULT_GLUT_EXCLUDE <- c("Glut-NSC", "Glut-NB", "Glut-Im", "Glut-GABA")

MIN_CELLS <- 10
LFC_MIN <- 0.25; PADJ_MAX <- 0.05; MIN_AUROC <- 0.70
TRANSGENE_TFS <- c("MAFB", "EMX2")   # never part of a signature: the OE target would score itself

# ---- pseudobulk ------------------------------------------------------------------------------------------------------
#' Sum counts per group. `f` is a per-cell grouping vector (NA cells dropped). Returns list(pb = genes x groups, n = cells per group).
pseudobulk_sum <- function(counts, f) {
  keep <- !is.na(f)
  f <- factor(f[keep]); counts <- counts[, keep, drop = FALSE]
  pb <- counts %*% sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  colnames(pb) <- levels(f)
  list(pb = pb, n = setNames(as.integer(table(f)), levels(f)))
}

#' DESeq2 on pseudobulk samples. `cd` is a data.frame (row names = pseudobulk ids) holding the design variables; `design`
#' a formula; the result is the `contrast`/`name` passed to results(). Rows kept when >= 5 counts in >= 2 samples, as in
#' figures/figure4/common.R. Returns gene, baseMean, lfc (unshrunk), lfc_se, stat, pvalue, padj.
deseq_pb <- function(pb, cd, design, contrast = NULL, name = NULL) {
  suppressMessages(library(DESeq2))
  mat <- round(as.matrix(pb[, rownames(cd), drop = FALSE]))
  dds <- DESeqDataSetFromMatrix(mat, cd, design)
  dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]
  dds <- DESeq(dds, quiet = TRUE)
  r <- if (!is.null(contrast)) results(dds, contrast = contrast) else results(dds, name = name)
  tibble(gene = rownames(r), baseMean = r$baseMean, lfc = r$log2FoldChange, lfc_se = r$lfcSE, stat = r$stat,
         pvalue = r$pvalue, padj = r$padj)
}

#' Use `full` when its model matrix is full rank on `cd`, else `reduced` (the figure4 rule for ~ lib + grp).
estimable_design <- function(cd, full, reduced) {
  mm <- model.matrix(full, cd)
  if (qr(mm)$rank == ncol(mm)) full else reduced
}

# ---- orthologs -------------------------------------------------------------------------------------------------------
#' OMA 1:1 chicken -> finch symbol map (the rule used by snrna/integration/: chicken gene_name_ck renamed to gene_name_zf).
oma_ck_to_zf <- function() {
  oma <- read_tsv(OMA_ORTHOLOGS_ZF_CK, show_col_types = FALSE, progress = FALSE) %>%
    filter(RelType == "1:1", !is.na(gene_name_ck), !is.na(gene_name_zf)) %>%
    distinct(gene_name_ck, gene_name_zf) %>%
    group_by(gene_name_ck) %>% filter(n() == 1) %>% group_by(gene_name_zf) %>% filter(n() == 1) %>% ungroup()
  setNames(oma$gene_name_zf, oma$gene_name_ck)
}

# ---- scoring helpers copied from songbird_dev nido/ (unchanged bodies) -----------------------------------------------
#' Mann-Whitney AUROC of `score` for `pos` vs the rest (calibrated_annotation_helpers.R).
fast_auroc <- function(score, pos) {
  stopifnot(length(score) == length(pos))
  pos <- as.logical(pos)
  keep <- !is.na(score) & !is.na(pos)
  score <- score[keep]; pos <- pos[keep]
  n_pos <- sum(pos); n_neg <- sum(!pos)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Equal-frequency expression bins, rank-based so all-zero ties still spread (calibrated_annotation_helpers.R).
bin_genes_by_expression <- function(stat, n_bins = 20) {
  stat <- stat[!is.na(stat)]
  br <- ceiling(rank(stat, ties.method = "first") / (length(stat) / n_bins))
  br[br > n_bins] <- n_bins
  setNames(as.integer(br), names(stat))
}

#' B random signatures matching the expression-bin composition (and UCell +/- signs) of `signature`
#' (calibrated_annotation_helpers.R).
expression_matched_sets <- function(signature, gene_bins, B = 100) {
  base <- sub("[-+]$", "", signature)
  signs <- ifelse(grepl("\\+$", signature), "+", ifelse(grepl("-$", signature), "-", ""))
  in_map <- base %in% names(gene_bins)
  base <- base[in_map]; signs <- signs[in_map]
  if (!length(base)) return(replicate(B, character(0), simplify = FALSE))
  want_bins <- gene_bins[base]
  by_bin <- split(names(gene_bins), gene_bins)
  lapply(seq_len(B), function(b) {
    drawn <- vapply(seq_along(want_bins), function(i) {
      pool <- by_bin[[as.character(want_bins[i])]]
      if (is.null(pool) || !length(pool)) return(NA_character_)
      pool[sample.int(length(pool), 1L)]
    }, character(1))
    ok <- !is.na(drawn) & !duplicated(drawn)
    paste0(drawn[ok], signs[ok])
  })
}

#' Seeded AUCell rankings: AUCell breaks ties at random and single-cell data is mostly ties at zero, so unseeded builds
#' differ run to run (core_signature_helpers.R).
build_aucell_rankings <- function(mat, seed = 10) {
  set.seed(seed)
  AUCell::AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE)
}

#' Up-only scoring of unsigned gene sets with one rule. `method` "ucell" (the songbird_dev default) or "aucell".
#' UCell maxRank 1500 (raised if a set is longer); AUCell aucMaxRank = 5% of features. Returns cells x sets.
score_sets <- function(mat, sets, method = c("ucell", "aucell"), rankings = NULL, ucell_maxrank = 1500, ncores = 1) {
  method <- match.arg(method)
  sets <- sets[lengths(sets) > 0]
  if (method == "ucell") {
    mr <- max(ucell_maxrank, max(lengths(sets)) + 1L)
    sc <- UCell::ScoreSignatures_UCell(mat, features = sets, maxRank = mr, ncores = ncores, name = "")
    sc[colnames(mat), names(sets), drop = FALSE]
  } else {
    if (is.null(rankings)) rankings <- build_aucell_rankings(mat)
    auc <- AUCell::AUCell_calcAUC(sets, rankings, aucMaxRank = ceiling(0.05 * nrow(mat)), verbose = FALSE)
    t(AUCell::getAUC(auc))[colnames(mat), names(sets), drop = FALSE]
  }
}
