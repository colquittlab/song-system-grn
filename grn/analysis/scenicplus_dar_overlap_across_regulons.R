#!/usr/bin/env Rscript
# Is the DAR overlap behind MAFB's RA / PVALB-1 genes unusual? The same statistic, for every other regulon.
#
# For MAFB the overlap of the RA-vs-C1H-1 and PVALB-1-vs-LAMP5 DAR sets among the regions behind its up-in-both genes (log2FC > 1 in both) was
# 1.2x what independence predicts, but independence is not a fair null: genes chosen for being up have more accessible regions whatever the
# TF. So the identical analysis is run for every TF with a +/+ regulon in config37 (all genes were tested in both contrasts):
#   whole regulon   regions linked to all of the TF's +/+ targets
#   up in both      regions linked to the TF's targets with log2FC > UP_LFC in both contrasts (the MAFB analysis)
# Statistics: fraction of the regions that are a DAR in both, observed / expected-if-independent, Jaccard of the two DAR sets, and a z-score
# against random same-size subsets of that TF's own targets (within-regulon null). MAFB is then ranked among the regulons.
#
#   Rscript scenicplus_dar_overlap_across_regulons.R
suppressMessages({library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
set.seed(2026)
MAIN <- "config37"
UP_LFC <- 1; MIN_REGIONS_ALL <- 100; MIN_GENES_UP <- 5; MIN_REGIONS_UP <- 30; NPERM <- 500
STORE <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/pycisTopic"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
bed_names <- function(p) { b <- read_tsv(p, col_names = c("c", "s", "e"), col_types = "cii", progress = FALSE); paste0(b$c, ":", b$s, "-", b$e) }
DAR_RA <- bed_names(file.path(STORE, "region_sets_k40_archr/DARs_song-pairs/Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.bed"))
DAR_PV <- bed_names(file.path(STORE, "dars_extra/GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.bed"))
cat("DARs: RA", length(DAR_RA), "| PVALB-1", length(DAR_PV), "\n")

e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(grepl("\\+/\\+$", eRegulon_name)) %>% distinct(TF, Gene, Region)
g <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE)
up_both <- g$gene[!is.na(g$lfc_RA_vs_C1H1) & !is.na(g$lfc_PV1_vs_low) & g$lfc_RA_vs_C1H1 > UP_LFC & g$lfc_PV1_vs_low > UP_LFC]
tested <- g$gene[!is.na(g$lfc_RA_vs_C1H1) & !is.na(g$lfc_PV1_vs_low)]
cat("genes tested in both contrasts:", length(tested), "| up in both (log2FC >", UP_LFC, "):", length(up_both), "\n")

stat <- function(regs) {
  n <- length(regs); a <- sum(regs %in% DAR_RA); b <- sum(regs %in% DAR_PV); both <- sum(regs %in% DAR_RA & regs %in% DAR_PV)
  ex <- a * b / n
  tibble(n_regions = n, n_RA = a, n_PV1 = b, n_both = both, frac_both = both / n, expected = ex, obs_over_exp = ifelse(ex > 0, both / ex, NA_real_),
         jaccard = ifelse(a + b - both > 0, both / (a + b - both), NA_real_))
}
res <- list()
for (tf in unique(e$TF)) {
  et <- e %>% filter(TF == tf, Gene %in% tested)
  if (n_distinct(et$Region) < MIN_REGIONS_ALL) next
  by_gene <- split(et$Region, et$Gene)
  all_regs <- unique(et$Region)
  res[[length(res) + 1]] <- bind_cols(TF = tf, analysis = "whole regulon", n_genes = length(by_gene), stat(all_regs), z_vs_random_targets = NA_real_)
  G <- intersect(names(by_gene), up_both)
  regs_G <- unique(unlist(by_gene[G]))
  if (length(G) >= MIN_GENES_UP && length(regs_G) >= MIN_REGIONS_UP) {
    obs <- stat(regs_G)
    nulls <- replicate(NPERM, { s <- sample(names(by_gene), length(G)); r <- unique(unlist(by_gene[s])); mean(r %in% DAR_RA & r %in% DAR_PV) })
    z <- (obs$frac_both - mean(nulls)) / sd(nulls)
    shared <- mean(sapply(by_gene[G], function(r) any(r %in% DAR_RA & r %in% DAR_PV)))
    res[[length(res) + 1]] <- bind_cols(TF = tf, analysis = "up in both", n_genes = length(G), obs, z_vs_random_targets = z) %>% mutate(frac_genes_sharing_a_DAR = shared)
  }
}
R <- bind_rows(res)
write_csv(R %>% mutate(across(where(is.double), ~ round(.x, 4))), file.path(HERE, paste0("scenicplus_dar_overlap_across_regulons_", MAIN, ".csv")))
genome <- length(intersect(DAR_RA, DAR_PV)) / (length(DAR_RA) * length(DAR_PV) / 499348)
cat(sprintf("genome-wide: %d regions are a DAR in both (%.2fx independence)\n", length(intersect(DAR_RA, DAR_PV)), genome))

for (an in c("whole regulon", "up in both")) {
  d <- R %>% filter(analysis == an)
  m <- d %>% filter(TF == "MAFB")
  cat(sprintf("\n== %s: %d regulons ==\n", an, nrow(d)))
  if (nrow(m)) {
    cat(sprintf("MAFB: %d genes, %d regions, %d in both DARs (frac %.3f), obs/exp %.2f, Jaccard %.3f%s\n", m$n_genes, m$n_regions, m$n_both, m$frac_both, m$obs_over_exp, m$jaccard,
                if (!is.na(m$z_vs_random_targets)) sprintf(", z vs random own targets %.1f", m$z_vs_random_targets) else ""))
    for (k in c("frac_both", "obs_over_exp", "jaccard")) cat(sprintf("  %-13s rank %d of %d (percentile %.0f); other regulons: median %.3f, IQR %.3f-%.3f, max %.3f\n", k, 1 + sum(d[[k]] > m[[k]], na.rm = TRUE), nrow(d),
                                                             100 * mean(d[[k]] <= m[[k]], na.rm = TRUE), median(d[[k]][d$TF != "MAFB"], na.rm = TRUE), quantile(d[[k]][d$TF != "MAFB"], .25, na.rm = TRUE),
                                                             quantile(d[[k]][d$TF != "MAFB"], .75, na.rm = TRUE), max(d[[k]][d$TF != "MAFB"], na.rm = TRUE)))
  } else cat("MAFB not in this analysis\n")
  print(as.data.frame(d %>% filter(TF != "MAFB") %>% arrange(desc(obs_over_exp)) %>% head(8) %>% select(TF, n_genes, n_regions, n_both, frac_both, obs_over_exp, jaccard, z_vs_random_targets) %>%
                        mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
}

# figure: distribution of obs/expected and of the fraction of regions that are DARs in both, across regulons, MAFB marked
P <- R %>% mutate(is_mafb = TF == "MAFB", analysis = factor(analysis, levels = c("whole regulon", "up in both")))
p <- ggplot(P, aes(analysis, obs_over_exp)) +
  geom_hline(yintercept = 1, linewidth = 0.3, color = FIG_INK_MUTED) +
  geom_jitter(data = P %>% filter(!is_mafb), width = 0.18, height = 0, size = 0.8, color = FIG_INK_MUTED, alpha = 0.7, stroke = 0) +
  geom_point(data = P %>% filter(is_mafb), size = 2, color = FIG_PAL[2], stroke = 0) +
  ggrepel::geom_text_repel(data = P %>% filter(is_mafb), aes(label = "MAFB"), nudge_x = 0.35, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                           segment.size = 0.2, min.segment.length = 0) +
  labs(x = NULL, y = "DAR overlap, observed / expected\n(RA vs C1H-1 and PVALB-1 vs LAMP5)") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank())
fig_save(p, file.path(OUT, MAIN, "dar_overlap_across_regulons"), width = 3, height = 2.3)
cat("saved figure\n")
