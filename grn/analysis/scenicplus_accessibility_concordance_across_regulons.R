#!/usr/bin/env Rscript
# Directional similarity of accessibility changes across regulons (no DAR calls; the PVALB-1 vs LAMP5 DAR test is underpowered).
#
# For each TF with a +/+ regulon in config37, take the regions linked to its targets and compare their ArchR log2FC in the two contrasts, using
# EVERY region whatever its significance:
#   rho            Spearman correlation of log2FC(RA vs C1H-1) and log2FC(PVALB-1 vs LAMP5) across the regulon's regions
#   both_up        fraction of regions with log2FC > 0 in both; ratio to the product of the two marginal fractions (what independence predicts)
#   concordance    fraction of regions where the two log2FCs have the same sign
# Two views, as for the DAR comparison (scenicplus_dar_overlap_across_regulons.R): the whole regulon, and the regions of the TF's targets that are up in
# both contrasts (RNA log2FC > UP_LFC in both; the MAFB analysis). For the second view a z-score compares rho with random same-size subsets of the
# TF's own targets. MAFB is ranked among the regulons; the genome-wide rho over all regions is the baseline.
#
#   Rscript scenicplus_accessibility_concordance_across_regulons.R
suppressMessages({library(SummarizedExperiment); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
set.seed(2026)
MAIN <- "config37"
UP_LFC <- 1; MIN_REGIONS_ALL <- 100; MIN_GENES_UP <- 5; MIN_REGIONS_UP <- 30; NPERM <- 500
STORE <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/archr_consensus"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")

get_lfc <- function(rds) {
  se <- readRDS(file.path(STORE, rds)); rr <- rowData(se)
  setNames(as.numeric(assay(se, "Log2FC")[, 1]), paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end))
}
LFC_RA <- get_lfc("markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds")
LFC_PV <- get_lfc("markers_GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.rds")
common <- intersect(names(LFC_RA), names(LFC_PV)); LFC_RA <- LFC_RA[common]; LFC_PV <- LFC_PV[common]
ok <- is.finite(LFC_RA) & is.finite(LFC_PV); LFC_RA <- LFC_RA[ok]; LFC_PV <- LFC_PV[ok]
rho_genome <- cor(LFC_RA, LFC_PV, method = "spearman")
cat(sprintf("regions with both log2FC: %d | genome-wide Spearman rho = %.3f | log2FC > 0 in both: %.3f (marginals %.3f x %.3f)\n", length(LFC_RA), rho_genome,
            mean(LFC_RA > 0 & LFC_PV > 0), mean(LFC_RA > 0), mean(LFC_PV > 0)))

e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(grepl("\\+/\\+$", eRegulon_name)) %>% distinct(TF, Gene, Region) %>% filter(Region %in% names(LFC_RA))
g <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE)
tested <- g$gene[!is.na(g$lfc_RA_vs_C1H1) & !is.na(g$lfc_PV1_vs_low)]
up_both <- g$gene[!is.na(g$lfc_RA_vs_C1H1) & !is.na(g$lfc_PV1_vs_low) & g$lfc_RA_vs_C1H1 > UP_LFC & g$lfc_PV1_vs_low > UP_LFC]

stat <- function(regs) {
  a <- LFC_RA[regs]; b <- LFC_PV[regs]
  pa <- mean(a > 0); pb <- mean(b > 0); both <- mean(a > 0 & b > 0)
  tibble(n_regions = length(regs), rho = suppressWarnings(cor(a, b, method = "spearman")), both_up = both, both_up_over_exp = both / (pa * pb),
         concordance = mean(sign(a) == sign(b)))
}
res <- list()
for (tf in unique(e$TF)) {
  et <- e %>% filter(TF == tf, Gene %in% tested)
  if (n_distinct(et$Region) < MIN_REGIONS_ALL) next
  by_gene <- split(et$Region, et$Gene)
  res[[length(res) + 1]] <- bind_cols(TF = tf, analysis = "whole regulon", n_genes = length(by_gene), stat(unique(et$Region)), z_rho_vs_random_targets = NA_real_)
  G <- intersect(names(by_gene), up_both); rG <- unique(unlist(by_gene[G]))
  if (length(G) >= MIN_GENES_UP && length(rG) >= MIN_REGIONS_UP) {
    obs <- stat(rG)
    nulls <- replicate(NPERM, { s <- sample(names(by_gene), length(G)); r <- unique(unlist(by_gene[s])); suppressWarnings(cor(LFC_RA[r], LFC_PV[r], method = "spearman")) })
    res[[length(res) + 1]] <- bind_cols(TF = tf, analysis = "up in both", n_genes = length(G), obs, z_rho_vs_random_targets = (obs$rho - mean(nulls, na.rm = TRUE)) / sd(nulls, na.rm = TRUE))
  }
}
R <- bind_rows(res)
write_csv(R %>% mutate(across(where(is.double), ~ round(.x, 4))), file.path(HERE, paste0("scenicplus_accessibility_concordance_across_regulons_", MAIN, ".csv")))

for (an in c("whole regulon", "up in both")) {
  d <- R %>% filter(analysis == an); m <- d %>% filter(TF == "MAFB"); o <- d %>% filter(TF != "MAFB")
  cat(sprintf("\n== %s: %d regulons ==\n", an, nrow(d)))
  cat(sprintf("MAFB: %d genes, %d regions | rho %.3f | both up %.3f (%.2fx marginals) | concordance %.3f%s\n", m$n_genes, m$n_regions, m$rho, m$both_up, m$both_up_over_exp, m$concordance,
              if (!is.na(m$z_rho_vs_random_targets)) sprintf(" | z vs random own targets %.1f", m$z_rho_vs_random_targets) else ""))
  for (k in c("rho", "both_up_over_exp", "concordance")) cat(sprintf("  %-17s MAFB %.3f | rank %d of %d | others: median %.3f, IQR %.3f to %.3f, max %.3f | regulons above MAFB: %.0f%%\n", k, m[[k]], 1 + sum(d[[k]] > m[[k]], na.rm = TRUE),
      nrow(d), median(o[[k]], na.rm = TRUE), quantile(o[[k]], .25, na.rm = TRUE), quantile(o[[k]], .75, na.rm = TRUE), max(o[[k]], na.rm = TRUE), 100 * mean(o[[k]] > m[[k]], na.rm = TRUE)))
  cat(sprintf("  regulons with rho > 0: %d of %d; with rho > genome-wide (%.3f): %d\n", sum(d$rho > 0, na.rm = TRUE), nrow(d), rho_genome, sum(d$rho > rho_genome, na.rm = TRUE)))
  print(as.data.frame(d %>% arrange(desc(rho)) %>% head(6) %>% select(TF, n_genes, n_regions, rho, both_up_over_exp, concordance, z_rho_vs_random_targets) %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
  cat("  most negative:", paste(d %>% arrange(rho) %>% head(4) %>% transmute(s = sprintf("%s (%.2f)", TF, rho)) %>% pull(s), collapse = ", "), "\n")
}

P <- R %>% mutate(is_mafb = TF == "MAFB", analysis = factor(analysis, levels = c("whole regulon", "up in both")))
p <- ggplot(P, aes(analysis, rho)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = FIG_INK_MUTED) +
  geom_hline(yintercept = rho_genome, linewidth = 0.3, linetype = "dashed", color = FIG_INK_MUTED) +
  geom_jitter(data = P %>% filter(!is_mafb), width = 0.18, height = 0, size = 0.8, color = FIG_INK_MUTED, alpha = 0.7, stroke = 0) +
  geom_point(data = P %>% filter(is_mafb), size = 2, color = FIG_PAL[2], stroke = 0) +
  ggrepel::geom_text_repel(data = P %>% filter(is_mafb), aes(label = "MAFB"), nudge_x = 0.35, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                           segment.size = 0.2, min.segment.length = 0) +
  labs(x = NULL, y = "Spearman rho of region log2FC\n(RA vs C1H-1, PVALB-1 vs LAMP5)") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank())
fig_save(p, file.path(OUT, MAIN, "accessibility_concordance_across_regulons"), width = 3, height = 2.3)
cat("\nsaved figure (dashed line: genome-wide rho)\n")
