#!/usr/bin/env Rscript
# Permutation test of the upper-right (UR) quadrant ratio per regulon (scenicplus_regulon_quadrant_fractions.R). Regulon labels are shuffled among the
# pooled genes (or regions) of all regulons: each permutation gives every regulon a random set of the same size drawn without replacement from the
# pool, and the UR ratio is recomputed for the random set (frac UR / (frac log2FC>0 in RA contrast x frac log2FC>0 in PVALB-1 contrast), both fractions
# taken from the set itself). This is a cross-regulon null: a regulon is unusual only relative to what arbitrary regulon-sized sets from the same pool
# give, so a regulon is not rewarded merely for being up in both contrasts for reasons shared by all targets. Reported per regulon: observed ratio and
# frac_UR, null mean, empirical one-sided p (upper tail; (1 + #null >= obs)/(1 + B)), z, BH q. Also the same for frac_UR itself.
# Genes belong to several regulons, so the pool is the union of regulon genes (regions: union of regulon regions); sets are not forced to be disjoint.
#
#   Rscript scenicplus_regulon_quadrant_permutation.R
suppressMessages({library(SummarizedExperiment); library(tidyverse); library(here)})
MAIN <- "config37"; MIN_REGIONS <- 100; B_GENE <- 10000; B_REGION <- 2000
ARCHR <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/archr_consensus"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
get_lfc <- function(rds) { se <- readRDS(file.path(ARCHR, rds)); rr <- rowData(se); setNames(as.numeric(assay(se, "Log2FC")[, 1]), paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end)) }
LFC_RA <- get_lfc("markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds"); LFC_PV <- get_lfc("markers_GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.rds")
g <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE) %>% filter(!is.na(lfc_RA_vs_C1H1), !is.na(lfc_PV1_vs_low))
e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(grepl("\\+/\\+$", eRegulon_name)) %>% distinct(TF, Gene, Region) %>% filter(Gene %in% g$gene, Region %in% names(LFC_RA), Region %in% names(LFC_PV))
gx <- setNames(g$lfc_RA_vs_C1H1, g$gene); gy <- setNames(g$lfc_PV1_vs_low, g$gene)

ratio_of <- function(a, b) { ur <- mean(a & b); ex <- mean(a) * mean(b); c(ur = ur, ratio = ifelse(ex > 0, ur / ex, NA_real_)) }
perm_kind <- function(sets, xa, yb, B, kind) {
  pool <- names(xa); A <- xa > 0; Bp <- yb > 0
  set.seed(20260)
  map_dfr(names(sets), function(tf) {
    idx <- match(sets[[tf]], pool); n <- length(idx)
    obs <- ratio_of(A[idx], Bp[idx])
    nul <- vapply(seq_len(B), function(i) { s <- sample.int(length(pool), n); ratio_of(A[s], Bp[s]) }, numeric(2))
    tibble(TF = tf, kind = kind, n = n, frac_UR = obs["ur"], ratio = obs["ratio"], null_ratio_mean = mean(nul["ratio", ], na.rm = TRUE), null_ratio_sd = sd(nul["ratio", ], na.rm = TRUE),
           null_frac_mean = mean(nul["ur", ]), null_frac_sd = sd(nul["ur", ]), p_ratio = (1 + sum(nul["ratio", ] >= obs["ratio"], na.rm = TRUE)) / (1 + B), p_frac = (1 + sum(nul["ur", ] >= obs["ur"])) / (1 + B))
  }) %>% mutate(z_ratio = (ratio - null_ratio_mean) / null_ratio_sd, z_frac = (frac_UR - null_frac_mean) / null_frac_sd, q_ratio = p.adjust(p_ratio, "BH"), q_frac = p.adjust(p_frac, "BH"),
                rank_ratio = rank(-ratio, ties.method = "min"), rank_frac = rank(-frac_UR, ties.method = "min"))
}
tfs <- e %>% group_by(TF) %>% filter(n_distinct(Region) >= MIN_REGIONS) %>% pull(TF) %>% unique() %>% sort()
gsets <- lapply(setNames(tfs, tfs), function(tf) unique(e$Gene[e$TF == tf])); rsets <- lapply(setNames(tfs, tfs), function(tf) unique(e$Region[e$TF == tf]))
gpool <- sort(unique(unlist(gsets))); rpool <- sort(unique(unlist(rsets)))
cat(sprintf("%d regulons; gene pool %d, region pool %d; B = %d / %d\n", length(tfs), length(gpool), length(rpool), B_GENE, B_REGION))
RG <- perm_kind(gsets, gx[gpool], gy[gpool], B_GENE, "gene_expression")
RR <- perm_kind(rsets, LFC_RA[rpool], LFC_PV[rpool], B_REGION, "region_accessibility")
R <- bind_rows(RG, RR) %>% mutate(across(where(is.double), ~ signif(.x, 4)))
write_csv(R, file.path(HERE, paste0("scenicplus_regulon_quadrant_permutation_", MAIN, ".csv")))
write_csv(R, file.path(OUT, MAIN, "regulon_scatter_plots", "quadrant_permutation.csv"))
for (k in unique(R$kind)) {
  d <- R %>% filter(kind == k); m <- d %>% filter(TF == "MAFB")
  cat(sprintf("\n== %s ==\nMAFB: n %d, frac_UR %.3f (null mean %.3f), ratio %.2f (null %.2f +/- %.2f, z %.2f), p_ratio %.4g, q %.3g, rank ratio %d/%d; p_frac %.4g, q %.3g, rank frac %d/%d\n", k, m$n, m$frac_UR,
              m$null_frac_mean, m$ratio, m$null_ratio_mean, m$null_ratio_sd, m$z_ratio, m$p_ratio, m$q_ratio, m$rank_ratio, nrow(d), m$p_frac, m$q_frac, m$rank_frac, nrow(d)))
  cat(sprintf("regulons with q_ratio < 0.05: %d (p < 0.05: %d of %d)\n", sum(d$q_ratio < 0.05), sum(d$p_ratio < 0.05), nrow(d)))
  cat("top by ratio z:", paste(d %>% arrange(desc(z_ratio)) %>% head(6) %>% transmute(s = sprintf("%s ratio %.2f z %.1f q %.2g n=%d", TF, ratio, z_ratio, q_ratio, n)) %>% pull(s), collapse = "; "), "\n")
}
