#!/usr/bin/env Rscript
# Do MAFB's target genes form a shared program in RA and PVALB-2 neurons?  (gene level; the region level is
# scenicplus_mafb_shared_regions.py, which reads the specificity table written here.)
#
# The notebook's scatter correlates two fold changes (RA vs CACNA1H-1, PVALB-2 vs LGE-1). That is tied to two arbitrary contrasts and
# its null (random genes) is wide (r = 0.11 +/- 0.13), so it cannot show reuse. Here each gene gets a specificity z-score for RA and
# for PVALB-2 across the other cell types, and the MAFB target set is compared with random gene sets MATCHED on mean expression and
# on how broadly the gene is expressed (so "expressed in all neurons" does not count as shared):
#   S1  number of targets with z > 1 in both        S2  mean z_RA * z_partner        S3  Spearman(z_RA, z_partner)
# z is computed across neuron clusters only (Glut-*/GABA-*), the stricter reference; the all-cluster version is also reported.
# The same statistic is computed for every other TF's +/+ regulon, to see whether MAFB stands out or any large regulon does.
#
#   Rscript scenicplus_mafb_shared_program.R
suppressMessages({library(Seurat); library(qs2); library(tidyverse); library(here)})
set.seed(2026)
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
OBJ <- "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
HERE <- here::here("grn/analysis")
RA <- "Glut-CACNA1H-RA"
## partner cell type (default PVALB-2); other partners write files with a suffix, e.g. Rscript ... GABA-MGE-SST-1
PV <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(PV)) PV <- "GABA-MGE-PVALB-2"
TAG <- gsub("-", "", sub("^GABA-MGE-", "", PV))
SUF <- if (PV == "GABA-MGE-PVALB-2") "" else paste0("_vs_", TAG)
NPERM <- 2000

obj <- qs_read(OBJ, nthreads = 8)
obj$cluster <- obj$cluster_hybrid
n_cells <- table(obj$cluster)
avg <- log1p(AverageExpression(obj, assays = "SCT", layer = "data", group.by = "cluster")[[1]])   # mean per cluster, size-independent
rm(obj); invisible(gc())
keep <- names(n_cells)[n_cells >= 50]
avg <- avg[, intersect(colnames(avg), keep)]
neu <- grep("^(Glut|GABA)-", colnames(avg), value = TRUE)
stopifnot(PV %in% colnames(avg), RA %in% colnames(avg))
cat("partner:", PV, "| clusters:", ncol(avg), "| neuron clusters:", length(neu), "\n")

zs <- function(m) { s <- apply(m, 1, sd); s[s == 0] <- NA; (m - rowMeans(m)) / s }
z_all <- zs(avg); z_neu <- zs(avg[, neu])
pool <- rownames(avg)[apply(avg, 1, max) >= 0.25 & !is.na(z_neu[, RA]) & !is.na(z_neu[, PV])]
mu <- rowMeans(avg[pool, neu]); sdv <- apply(avg[pool, neu], 1, sd)
strata <- paste(cut(mu, quantile(mu, 0:10 / 10), include.lowest = TRUE, labels = FALSE),
                cut(sdv, quantile(sdv, 0:4 / 4), include.lowest = TRUE, labels = FALSE), sep = "_")
names(strata) <- pool
cat("expressed pool:", length(pool), "genes in", length(unique(strata)), "expression/breadth strata\n")
write_csv(tibble(gene = pool, mean_expr_neurons = mu, sd_expr_neurons = sdv, stratum = strata,
                 z_RA = z_neu[pool, RA], z_partner = z_neu[pool, PV], zall_RA = z_all[pool, RA], zall_partner = z_all[pool, PV],
                 expr_RA = avg[pool, RA], expr_partner = avg[pool, PV]), file.path(OUT, paste0("mafb_gene_specificity_RA_", if (SUF == "") "PVALB2" else TAG, ".csv")))

stat <- function(g, z) { a <- z[g, RA]; b <- z[g, PV]; c(S1 = sum(a > 1 & b > 1), S2 = mean(a * b), S3 = suppressWarnings(cor(a, b, method = "spearman"))) }
enrich <- function(targets, z = z_neu, nperm = NPERM) {
  g <- intersect(targets, pool)
  if (length(g) < 20) return(NULL)
  obs <- stat(g, z)
  by_s <- split(setdiff(pool, g), strata[setdiff(pool, g)])
  need <- table(strata[g])
  null <- replicate(nperm, {
    s <- unlist(lapply(names(need), function(k) { p <- by_s[[k]]; if (is.null(p) || !length(p)) character() else sample(p, need[[k]], replace = length(p) < need[[k]]) }))
    stat(s, z)
  })
  tibble(n_targets = length(g), S1 = obs["S1"], S1_null = mean(null["S1", ]), S1_ratio = obs["S1"] / mean(null["S1", ]),
         S1_p = (1 + sum(null["S1", ] >= obs["S1"])) / (nperm + 1),
         S2 = obs["S2"], S2_null = mean(null["S2", ]), S2_z = (obs["S2"] - mean(null["S2", ])) / sd(null["S2", ]),
         S3 = obs["S3"], S3_null = mean(null["S3", ]), S3_z = (obs["S3"] - mean(null["S3", ])) / sd(null["S3", ]))
}
ldc <- function(cfg) {
  p <- file.path(OUT, cfg, "scenicplus_eRegulons.txt")
  if (!file.exists(p)) return(NULL)
  read_tsv(p, col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c"), progress = FALSE) %>% filter(grepl("\\+/\\+$", eRegulon_name))
}

# 1. MAFB +/+ targets, per config (the correct-input configs of the sweep, plus the 40-topic runs)
cfgs <- paste0("config", c(1:11, 16:28, 34:39))
res <- list(); mafb <- list()
for (cfg in cfgs) {
  d <- ldc(cfg); if (is.null(d)) next
  g <- unique(d$Gene[d$TF == "MAFB"]); mafb[[cfg]] <- g
  r <- enrich(g); if (!is.null(r)) res[[cfg]] <- bind_cols(set = cfg, r)
}
# 2. targets present in at least half of the sweep configs (the stable core)
mem <- read_csv(file.path(HERE, "scenicplus_mafb_target_membership.csv"), show_col_types = FALSE)
core <- mem$Gene[mem$n_pp >= 12]
res[["stable core (>=12 of 24 configs)"]] <- bind_cols(set = "stable core (>=12 of 24 configs)", enrich(core))
# 3. same, all-cluster z (broader reference)
res[["config1, z across all clusters"]] <- bind_cols(set = "config1, z across all clusters", enrich(mafb[["config1"]], z = z_all))
tab <- bind_rows(res)
write_csv(tab, file.path(HERE, paste0("scenicplus_mafb_shared_program", SUF, ".csv")))
print(as.data.frame(tab %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# 4. every other TF's +/+ regulon in config1: does MAFB stand out?
d1 <- ldc("config1")
sets <- d1 %>% group_by(TF) %>% summarise(genes = list(unique(Gene)), n = n_distinct(Gene)) %>% filter(n >= 50)
ref <- map2_dfr(sets$TF, sets$genes, function(tf, g) { r <- enrich(g, nperm = 500); if (is.null(r)) NULL else bind_cols(TF = tf, r) })
write_csv(ref, file.path(HERE, paste0("scenicplus_shared_program_all_tfs_config1", SUF, ".csv")))
cat("\nconfig1: MAFB vs", nrow(ref) - 1, "other TFs with >= 50 +/+ targets\n")
for (m in c("S1_ratio", "S2_z", "S3_z")) cat(sprintf("  %s: MAFB %.2f, rank %d of %d (higher = more shared); other TFs median %.2f (IQR %.2f to %.2f)\n", m,
  ref[[m]][ref$TF == "MAFB"], rank(-ref[[m]])[ref$TF == "MAFB"], nrow(ref), median(ref[[m]][ref$TF != "MAFB"]),
  quantile(ref[[m]][ref$TF != "MAFB"], .25), quantile(ref[[m]][ref$TF != "MAFB"], .75)))
print(as.data.frame(ref %>% arrange(desc(S2_z)) %>% head(8) %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% select(TF, n_targets, S1, S1_ratio, S2_z, S3_z)), row.names = FALSE)
