#!/usr/bin/env Rscript
# Share of each regulon's genes (RNA log2FC) and regions (accessibility log2FC) in the upper-right quadrant of the RA vs C1H-1 / PVALB-1 vs LAMP5
# scatters (scenicplus_regulon_scatter_plots.R): log2FC > 0 in both contrasts. Reported per regulon and per version:
#   frac_UR           n in the upper-right quadrant / n
#   frac_expected     product of the two marginal fractions (log2FC > 0 in RA vs C1H-1) x (log2FC > 0 in PVALB-1 vs LAMP5), i.e. what independence predicts
#   ratio             frac_UR / frac_expected  (a regulon can have a large frac_UR simply because it is up in both contrasts; the ratio corrects for that)
#   fisher_p          Fisher exact test of the 2x2 table (sign in RA contrast x sign in PVALB-1 contrast); odds ratio alongside
#   frac_UR_strict    share with log2FC > STRICT in both
# Same regulons and points as the scatter plots.
#
#   Rscript scenicplus_regulon_quadrant_fractions.R
suppressMessages({library(SummarizedExperiment); library(tidyverse); library(here)})
MAIN <- "config37"; MIN_REGIONS <- 100; STRICT <- 1
ARCHR <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/archr_consensus"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
get_lfc <- function(rds) { se <- readRDS(file.path(ARCHR, rds)); rr <- rowData(se); setNames(as.numeric(assay(se, "Log2FC")[, 1]), paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end)) }
LFC_RA <- get_lfc("markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds"); LFC_PV <- get_lfc("markers_GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.rds")
g <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE) %>% filter(!is.na(lfc_RA_vs_C1H1), !is.na(lfc_PV1_vs_low))
e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(grepl("\\+/\\+$", eRegulon_name)) %>% distinct(TF, Gene, Region) %>% filter(Gene %in% g$gene, Region %in% names(LFC_RA), Region %in% names(LFC_PV))
gx <- setNames(g$lfc_RA_vs_C1H1, g$gene); gy <- setNames(g$lfc_PV1_vs_low, g$gene)

quad <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  tab <- table(factor(x > 0, c(TRUE, FALSE)), factor(y > 0, c(TRUE, FALSE)))
  f <- fisher.test(tab)
  ur <- mean(x > 0 & y > 0); ex <- mean(x > 0) * mean(y > 0)
  tibble(n = n, n_UR = sum(x > 0 & y > 0), frac_UR = ur, frac_expected = ex, ratio = ur / ex, odds_ratio = unname(f$estimate), fisher_p = f$p.value,
         frac_UL = mean(x <= 0 & y > 0), frac_LR = mean(x > 0 & y <= 0), frac_LL = mean(x <= 0 & y <= 0), frac_UR_strict = mean(x > STRICT & y > STRICT))
}
res <- list()
for (tf in sort(unique(e$TF))) {
  et <- e %>% filter(TF == tf)
  if (n_distinct(et$Region) < MIN_REGIONS) next
  ge <- unique(et$Gene); rg <- unique(et$Region)
  res[[length(res) + 1]] <- bind_cols(TF = tf, kind = "gene_expression", quad(gx[ge], gy[ge]))
  res[[length(res) + 1]] <- bind_cols(TF = tf, kind = "region_accessibility", quad(LFC_RA[rg], LFC_PV[rg]))
}
R <- bind_rows(res) %>% group_by(kind) %>% mutate(fisher_q = p.adjust(fisher_p, "BH")) %>% ungroup()
write_csv(R %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(HERE, paste0("scenicplus_regulon_quadrant_fractions_", MAIN, ".csv")))
write_csv(R %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(OUT, MAIN, "regulon_scatter_plots", "quadrant_fractions.csv"))

for (k in unique(R$kind)) {
  d <- R %>% filter(kind == k); m <- d %>% filter(TF == "MAFB"); o <- d %>% filter(TF != "MAFB")
  cat(sprintf("\n== %s: %d regulons ==\n", k, nrow(d)))
  cat(sprintf("all regulons: frac_UR median %.3f (IQR %.3f-%.3f, range %.3f-%.3f) | ratio median %.2f (IQR %.2f-%.2f) | Fisher q < 0.05: %d positive, %d negative\n", median(d$frac_UR),
              quantile(d$frac_UR, .25), quantile(d$frac_UR, .75), min(d$frac_UR), max(d$frac_UR), median(d$ratio), quantile(d$ratio, .25), quantile(d$ratio, .75),
              sum(d$fisher_q < 0.05 & d$odds_ratio > 1), sum(d$fisher_q < 0.05 & d$odds_ratio < 1)))
  cat(sprintf("MAFB: n %d, upper-right %d (%.3f; expected %.3f; ratio %.2f; odds ratio %.2f, Fisher p %.3g) | frac_UR rank %d of %d, ratio rank %d of %d | strict (>%g in both) %.3f\n",
              m$n, m$n_UR, m$frac_UR, m$frac_expected, m$ratio, m$odds_ratio, m$fisher_p, 1 + sum(d$frac_UR > m$frac_UR), nrow(d), 1 + sum(d$ratio > m$ratio), nrow(d), STRICT, m$frac_UR_strict))
  cat("highest frac_UR:", paste(d %>% arrange(desc(frac_UR)) %>% head(6) %>% transmute(s = sprintf("%s %.2f (ratio %.2f, n=%d)", TF, frac_UR, ratio, n)) %>% pull(s), collapse = "; "), "\n")
  cat("highest ratio:  ", paste(d %>% arrange(desc(ratio)) %>% head(6) %>% transmute(s = sprintf("%s %.2f (frac %.2f, q=%.2g, n=%d)", TF, ratio, frac_UR, fisher_q, n)) %>% pull(s), collapse = "; "), "\n")
  cat("lowest ratio:   ", paste(d %>% arrange(ratio) %>% head(4) %>% transmute(s = sprintf("%s %.2f (frac %.2f, n=%d)", TF, ratio, frac_UR, n)) %>% pull(s), collapse = "; "), "\n")
}
