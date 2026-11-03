#!/usr/bin/env Rscript
# 04_signature_scoring.R -- do plasmid+ chicken glutamatergic neurons score more song-neuron-like than their neighbors?
#
# The scoring rule is the one settled on in songbird_dev (session "glutamatergic neuron clustering seurat", nido/):
# UCell, up-only sets, expression-matched permutation null (20 expression bins, B = 1000 sets per signature, scored with the
# identical rule); AUCell (seeded rankings, aucMaxRank = 5% of features) is run as the second scorer. Sets: for each axis
# (RA vs CACNA1H-1, HVCra vs DACH2-1, HVCx vs DACH2-4), SONG (genes up in the song neuron) and SISTER (genes up in its
# non-song sister), STRICT (adult-AUROC gated) and ALL. An RA-like shift is SONG up and SISTER down.
#
# What is NOT carried over: the adult-anchor calibration (0 = sister median, 1 = song median). Anchors measured on adult
# finch 10x nuclei do not transfer to embryonic chicken PIPseq nuclei (different species, platform, depth, gene
# detection), so the claim here is relative -- plasmid+ vs plasmid- cells of the same embryo -- which is what the
# songbird_dev doc recommends for an OE design.
#
# Inference, embryo as the unit:
#   - per-cell scores are adjusted for depth (residuals of score ~ library + log10 UMI over all glutamatergic cells), since
#     plasmid detection selects deeper cells;
#   - per embryo with >= MIN_POS plasmid+ glutamatergic cells: shift = mean adjusted score (plasmid+) - mean (plasmid-);
#   - per library: mean shift over its embryos; MAFB difference-in-differences: +dox minus -dox;
#   - each statistic is compared with the same statistic for the 1000 expression-matched null sets of that signature:
#     z = (obs - mean(null)) / sd(null), two-sided empirical p. This removes any shift that random genes of the same
#     expression level also show (depth, cell size, maturation of the electroporated cohort).
#
#   Rscript overexpression/reanalysis/04_signature_scoring.R
suppressMessages({library(Seurat); library(qs2)})
source(here::here("overexpression/reanalysis/oe_helpers.R"))
stopifnot(requireNamespace("UCell", quietly = TRUE), requireNamespace("AUCell", quietly = TRUE))
set.seed(2026)
MIN_POS <- 3; B <- 1000; N_BINS <- 20; NCORES <- 16
# tg+ definition (OE_TG_DEF), as in 03_pseudobulk_de.R: "plasmid" = any plasmid read (default); "inducible" = any MAFB-tg /
# EMX2-tg / mScarlet read, tg- = no plasmid read, and cells whose only plasmid read is the constitutive mNG are left out of
# the embryo means (they are still scored and still inform the depth model). Outputs of "inducible" carry _inducible.
TG_DEF <- match.arg(Sys.getenv("OE_TG_DEF", "plasmid"), c("plasmid", "inducible"))
SFX <- if (TG_DEF == "inducible") "_inducible" else ""
cat("tg+ definition:", TG_DEF, "\n")
MIN_SET <- 10   # sets with fewer genes present in the chicken matrix are dropped

obj <- qs_read(file.path(OE_OUT, "oe_obj.qs2"), nthreads = 8)
md <- obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>% filter(cell_class == "Glut")
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")[, md$cell]
rm(obj); invisible(gc())

## chicken -> finch gene names (OMA 1:1), log-normalized
ck2zf <- oma_ck_to_zf()
cnt <- cnt[intersect(rownames(cnt), names(ck2zf)), ]
rownames(cnt) <- ck2zf[rownames(cnt)]
lc <- log1p(t(t(cnt) / md$nCount_RNA) * 1e4)
cat(sprintf("chicken glutamatergic cells: %d; genes with a 1:1 finch ortholog: %d\n", ncol(lc), nrow(lc)))

## signatures
sig <- read_csv(file.path(OE_TABLES, "reference_signature_genes.csv"), show_col_types = FALSE) %>% filter(has_ortholog, !excluded_tf)
sets <- c(sig %>% filter(specific) %>% split(paste0(.$axis, "__", .$direction, "_STRICT")) %>% map("gene"),
          sig %>% split(paste0(.$axis, "__", .$direction, "_ALL")) %>% map("gene"))
sets <- map(sets, ~ intersect(.x, rownames(lc)))
detect <- rowMeans(lc > 0)
set_info <- tibble(set = names(sets), n_genes = lengths(sets), n_detected_1pct = map_int(sets, ~ sum(detect[.x] >= 0.01)))
cat("\nsignature sets in the chicken matrix:\n"); print(as.data.frame(set_info), row.names = FALSE)
sets <- sets[lengths(sets) >= MIN_SET]

## expression-matched null sets
gene_bins <- bin_genes_by_expression(rowMeans(lc), n_bins = N_BINS)
null_sets <- imap(sets, function(s, nm) setNames(expression_matched_sets(s, gene_bins, B = B), sprintf("%s__null%04d", nm, seq_len(B))))
all_sets <- c(sets, unlist(unname(null_sets), recursive = FALSE))
cat(sprintf("\nscoring %d real + %d null sets\n", length(sets), length(all_sets) - length(sets)))

## per-embryo statistics
md <- md %>% mutate(pos = if (TG_DEF == "inducible") inducible else plasmid,
                    use = if (TG_DEF == "inducible") inducible | !plasmid else TRUE,
                    logumi = log10(nCount_RNA))
X <- model.matrix(~ lib + logumi, md)
emb_ok <- md %>% filter(use) %>% group_by(embryo, lib, dox) %>% summarise(n_pos = sum(pos), n_neg = sum(!pos), .groups = "drop") %>% filter(n_pos >= MIN_POS, n_neg >= MIN_CELLS)
cat("embryos used:", paste0(emb_ok$embryo, " (", emb_ok$n_pos, "+/", emb_ok$n_neg, "-)", collapse = ", "), "\n")

stats_from_scores <- function(S) {                    # S: cells x sets
  R <- qr.resid(qr(X), S)                             # depth/library-adjusted scores
  shift <- sapply(emb_ok$embryo, function(e) {
    i <- md$embryo == e & md$use
    colMeans(R[i & md$pos, , drop = FALSE]) - colMeans(R[i & !md$pos, , drop = FALSE])
  })                                                  # sets x embryos
  if (is.null(dim(shift))) shift <- matrix(shift, nrow = ncol(S), dimnames = list(colnames(S), emb_ok$embryo))
  lib_mean <- sapply(split(emb_ok$embryo, emb_ok$lib), function(es) rowMeans(shift[, es, drop = FALSE]))
  out <- cbind(shift, lib_mean)
  if (all(c("MAFBp", "MAFBn") %in% colnames(lib_mean))) out <- cbind(out, MAFB_dox_DiD = lib_mean[, "MAFBp"] - lib_mean[, "MAFBn"])
  out
}

run_method <- function(method) {
  S <- if (method == "ucell") score_sets(lc, all_sets, "ucell", ncores = NCORES) else score_sets(lc, all_sets, "aucell", rankings = build_aucell_rankings(lc))
  st <- stats_from_scores(S)
  map_dfr(names(sets), function(nm) {
    nul <- st[grep(paste0("^", nm, "__null"), rownames(st)), , drop = FALSE]
    obs <- st[nm, ]
    mu <- colMeans(nul); sdv <- apply(nul, 2, sd)
    p2 <- sapply(seq_along(obs), function(j) (1 + sum(abs(nul[, j] - mu[j]) >= abs(obs[j] - mu[j]))) / (nrow(nul) + 1))
    tibble(method = method, set = nm, statistic = names(obs), obs = obs, null_mean = mu, null_sd = sdv, z = (obs - mu) / sdv, p_emp = p2)
  })
}
res <- bind_rows(run_method("ucell"), run_method("aucell")) %>%
  separate(set, c("axis", "set"), sep = "__", remove = FALSE) %>%
  mutate(level = case_when(statistic %in% emb_ok$embryo ~ "embryo", statistic == "MAFB_dox_DiD" ~ "DiD", TRUE ~ "library"))
write_csv(res %>% mutate(across(c(obs, null_mean, null_sd), ~ signif(.x, 4)), z = round(z, 3), p_emp = signif(p_emp, 3)),
          file.path(OE_TABLES, paste0("oe_signature_scores", SFX, ".csv")))
write_csv(set_info, file.path(OE_TABLES, "oe_signature_sets.csv"))

cat("\nlibrary-level and DiD shifts vs the expression-matched null (STRICT sets; z > 0 = higher in tg+ cells; tg+ = ", TG_DEF, "):\n", sep = "")
print(as.data.frame(res %>% filter(level != "embryo", grepl("STRICT", set)) %>% select(method, axis, set, statistic, z, p_emp) %>%
                      mutate(z = round(z, 2), p_emp = signif(p_emp, 2)) %>% pivot_wider(names_from = method, values_from = c(z, p_emp))), row.names = FALSE)
cat("\nper-embryo z (UCell, STRICT):\n")
print(as.data.frame(res %>% filter(level == "embryo", method == "ucell", grepl("STRICT", set)) %>% select(axis, set, statistic, z) %>%
                      mutate(z = round(z, 2)) %>% pivot_wider(names_from = statistic, values_from = z)), row.names = FALSE)

## figure: per-embryo z by library, faceted by axis; SONG and SISTER sets
pd <- res %>% filter(method == "ucell", level == "embryo", grepl("STRICT", set)) %>%
  left_join(emb_ok %>% select(statistic = embryo, lib, n_pos), by = "statistic") %>%
  mutate(lib = factor(lib, levels = c("MAFBn", "MAFBp", "EMX2p"), labels = c("MAFB -dox", "MAFB +dox", "EMX2 +dox")),
         axis = factor(axis, levels = AXES$axis), set = factor(sub("_STRICT", "", set), levels = c("SONG", "SISTER"), labels = c("song-up genes", "sister-up genes")))
pl <- ggplot(pd, aes(lib, z, color = set)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = FIG_GRID, alpha = 0.6) +
  geom_hline(yintercept = 0, color = FIG_AXIS, linewidth = 0.3) +
  geom_point(aes(size = n_pos), position = position_dodge(width = 0.5), alpha = 0.85) +
  scale_size_area(max_size = 2.2, breaks = c(5, 10, 18), name = if (TG_DEF == "inducible") "inducible+ cells" else "plasmid+ cells") +
  scale_color_manual(values = c(FIG_PAL[2], FIG_PAL[1]), name = NULL) +
  facet_wrap(~ axis, nrow = 1) + scale_y_continuous(expand = fig_expand()) +
  labs(x = NULL, y = "shift vs matched null (z)") + theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = FIG_PT_AXIS_TITLE), legend.position = "bottom",
        legend.text = element_text(size = FIG_PT_AXIS_TEXT), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.key.size = unit(6, "pt"),
        panel.grid.major.x = element_blank())
fig_save(pl, file.path(OE_OUT, paste0("signature_shift_per_embryo", SFX)), width = 4, height = 2.4)
cat("\ndone\n")
