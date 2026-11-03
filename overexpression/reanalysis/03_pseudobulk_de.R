#!/usr/bin/env Rscript
# 03_pseudobulk_de.R -- transgene detection, pseudobulk DE in glutamatergic neurons, and song-axis enrichment.
#
# 1. Detection. How often does a cell that carries/expresses the plasmid show zero plasmid reads (so "tg-" cells include
#    true positives, the reviewer's point), and how many "tg+" calls could be ambient RNA? A depth-offset zero-truncated
#    Poisson is fitted to the plasmid UMI counts of positive cells (rate per UMI of library depth); P(miss) for a cell of
#    depth d is exp(-rate * d). Ambient: souporcell's ambient fraction x the library's plasmid share of UMIs gives each
#    cell's expected ambient plasmid UMIs, hence the expected number of ambient-only positives.
#    Hidden positives in the tg- group ATTENUATE a tg+ vs tg- difference (they dilute it toward zero); they do not create
#    one. What can create one is that detection selects deeper cells (02_process_oe.R), so every within-library contrast
#    is paired with the same contrast in the -dox library, where the selection is identical but MAFB is not induced.
#
# 2. Pseudobulk DESeq2, glutamatergic cells, replicate = embryo (library x souporcell cluster). Plasmid+ = any plasmid read.
#    The project floor is >= 10 cells per pseudobulk sample; only one embryo has 10 plasmid+ glutamatergic cells, so the
#    plasmid+ samples here use >= MIN_POS (3) cells and the plasmid- samples >= 10. This is a deliberate, low-power
#    exception, written into every output.
#      within_<lib>    plasmid+ vs plasmid-, paired within embryo:   ~ embryo + grp
#      mafb_dox_did    (plasmid+ vs -) in +dox minus the same in -dox: ~ dox + dox:emb + dox:grp  (DESeq2 nested design)
#      mafb_dox_all    all glutamatergic cells, +dox vs -dox embryos: ~ dox
#
# 3. fgsea of each adult song-axis signature (01_reference_signatures.R) on each DE ranking (DESeq2 Wald stat), chicken
#    genes renamed to finch by OMA 1:1. RA-like = SONG sets up (NES > 0) and SISTER sets down (NES < 0).
#
#   Rscript overexpression/reanalysis/03_pseudobulk_de.R
suppressMessages({library(Seurat); library(qs2); library(DESeq2); library(fgsea)})
source(here::here("overexpression/reanalysis/oe_helpers.R"))
set.seed(2026)
MIN_POS <- 3
# tg+ definition (OE_TG_DEF): "plasmid" = any plasmid read (default); "inducible" = any dox-inducible read (MAFB-tg / EMX2-tg
# or mScarlet), with tg- = no plasmid read of any kind (cells whose only plasmid read is the constitutive mNG are dropped from
# both groups). Outputs of "inducible" carry _inducible.
TG_DEF <- match.arg(Sys.getenv("OE_TG_DEF", "plasmid"), c("plasmid", "inducible"))
SFX <- if (TG_DEF == "inducible") "_inducible" else ""
cat("tg+ definition:", TG_DEF, "\n")

obj <- qs_read(file.path(OE_OUT, "oe_obj.qs2"), nthreads = 8)
md <- obj@meta.data %>% rownames_to_column("cell") %>% as_tibble()
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
rm(obj); invisible(gc())

## 1. detection -------------------------------------------------------------------------------------------------------
ztp_rate <- function(y, d) {          # MLE of rate r with lambda_i = r * d_i, y_i >= 1 (zero-truncated Poisson)
  nll <- function(lr) { l <- exp(lr) * d; -sum(y * log(l) - l - log1p(-exp(-l))) }
  exp(optimize(nll, c(log(1e-8), log(1)))$minimum)
}
amb <- map_dbl(c(MAFBn = "bGalGal1.mat.mafb_v2/MAFBn/souporcell_2", MAFBp = "bGalGal1.mat.mafb_v2/MAFBp/souporcell_3",
                 EMX2p = "bGalGal1.mat.emx2/EMX2p/souporcell_4"),
               ~ as.numeric(str_extract(read_lines(file.path(PIPSEQ_MAFB_EMX2_DIR, .x, "ambient_rna.txt"))[1], "[0-9.]+(?=%)")) / 100)  # "ambient RNA estimated as 19.85%"
cat("souporcell ambient RNA fraction:", paste(names(amb), signif(amb, 3), collapse = ", "), "\n")

md <- md %>% mutate(pl_umi = mNG + mScarlet + tf_tg)
det <- md %>% group_by(lib, dox, tf) %>% group_modify(function(d, k) {
  pos <- d %>% filter(pl_umi > 0)
  r <- ztp_rate(pos$pl_umi, pos$nCount_RNA)
  glut <- d %>% filter(cell_class == "Glut")
  share <- sum(d$pl_umi) / sum(d$nCount_RNA)                       # plasmid share of this library's UMIs
  amb_pos <- sum(1 - exp(-amb[[k$lib]] * d$nCount_RNA * share))   # expected cells positive from ambient alone
  p_miss_glut_neg <- exp(-r * glut$nCount_RNA[!glut$plasmid])
  obs_frac <- mean(glut$plasmid)
  true_frac <- min(1, obs_frac / (1 - mean(exp(-r * glut$nCount_RNA))))
  tibble(cells = nrow(d), plasmid_pos = nrow(pos), mean_umi_in_pos = mean(pos$pl_umi), rate_per_10k_umi = 1e4 * r,
         p_miss_median_glut = median(p_miss_glut_neg), expected_ambient_pos = amb_pos, ambient_pos_frac = amb_pos / nrow(pos),
         glut_obs_pos_frac = obs_frac, glut_implied_true_frac = true_frac,
         glut_neg_contamination = (true_frac - obs_frac) / (1 - obs_frac))
}) %>% ungroup() %>% mutate(across(where(is.double), ~ signif(.x, 3)))
cat("\nplasmid detection (rate fitted on positive cells; P(miss) for plasmid- glutamatergic cells at their depth):\n")
print(as.data.frame(det), row.names = FALSE)
if (TG_DEF == "plasmid") write_csv(det, file.path(OE_TABLES, "oe_plasmid_detection.csv"))

## 2. pseudobulk ------------------------------------------------------------------------------------------------------
g <- md %>% filter(cell_class == "Glut")
g <- if (TG_DEF == "inducible") g %>% filter(inducible | !plasmid) %>% mutate(grp = ifelse(inducible, "pos", "neg")) else g %>% mutate(grp = ifelse(plasmid, "pos", "neg"))
pbs <- pseudobulk_sum(cnt[, g$cell], paste(g$embryo, g$grp, sep = "|"))
info <- tibble(id = names(pbs$n), n = pbs$n) %>% separate(id, c("embryo", "grp"), sep = "\\|", remove = FALSE) %>%
  left_join(distinct(g, embryo, lib, dox), by = "embryo") %>%
  mutate(ok = ifelse(grp == "pos", n >= MIN_POS, n >= MIN_CELLS))
# an embryo enters a paired contrast only when both of its samples pass
paired <- info %>% group_by(embryo) %>% filter(all(ok), n() == 2) %>% ungroup()
cat("\npseudobulk samples (glutamatergic; plasmid+ >= ", MIN_POS, " cells, plasmid- >= ", MIN_CELLS, "):\n", sep = "")
print(as.data.frame(info %>% select(embryo, grp, n, ok) %>% pivot_wider(names_from = grp, values_from = c(n, ok))), row.names = FALSE)
write_csv(info %>% select(-id), file.path(OE_TABLES, paste0("oe_pseudobulk_samples", SFX, ".csv")))

mk_cd <- function(s) data.frame(row.names = s$id, embryo = factor(s$embryo), lib = factor(s$lib), dox = factor(s$dox, levels = c("minus", "plus")),
                                grp = factor(s$grp, levels = c("neg", "pos")))
res <- list()
for (L in unique(paired$lib)) {
  s <- paired %>% filter(lib == L)
  if (n_distinct(s$embryo) < 2) { cat("  within_", L, ": < 2 paired embryos, skipped\n", sep = ""); next }
  res[[paste0("within_", L)]] <- deseq_pb(pbs$pb, mk_cd(s), ~ embryo + grp, contrast = c("grp", "pos", "neg"))
  cat(sprintf("  within_%s: %d embryos paired, %d genes\n", L, n_distinct(s$embryo), nrow(res[[paste0("within_", L)]])))
}
# difference-in-differences: nested-embryo design from the DESeq2 vignette ("group-specific condition effects, individuals nested within groups")
s <- paired %>% filter(lib %in% c("MAFBn", "MAFBp")) %>% group_by(dox) %>% mutate(emb = factor(as.integer(factor(embryo)))) %>% ungroup()
cd <- mk_cd(s); cd$emb <- s$emb
mm <- model.matrix(~ dox + dox:emb + dox:grp, cd); mm <- mm[, colSums(mm != 0) > 0]
dds <- DESeqDataSetFromMatrix(round(as.matrix(pbs$pb[, rownames(cd)])), cd, ~ 1)
dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]
dds <- DESeq(dds, full = mm, quiet = TRUE)
rn <- resultsNames(dds); cat("  mafb_dox_did coefficients:", paste(rn, collapse = ", "), "\n")
tidy_r <- function(r) tibble(gene = rownames(r), baseMean = r$baseMean, lfc = r$log2FoldChange, lfc_se = r$lfcSE, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
res[["mafb_dox_did"]] <- tidy_r(results(dds, contrast = list("doxplus.grppos", "doxminus.grppos")))
res[["mafb_plus_within_nested"]] <- tidy_r(results(dds, name = "doxplus.grppos"))
res[["mafb_minus_within_nested"]] <- tidy_r(results(dds, name = "doxminus.grppos"))
# whole-library dox effect on all glutamatergic cells
pb_all <- pseudobulk_sum(cnt[, g$cell], g$embryo)
s <- tibble(id = names(pb_all$n), embryo = id, n = pb_all$n) %>% left_join(distinct(g, embryo, lib, dox), by = "embryo") %>%
  filter(lib %in% c("MAFBn", "MAFBp"), n >= MIN_CELLS) %>% mutate(grp = "neg")
res[["mafb_dox_all"]] <- deseq_pb(pb_all$pb, mk_cd(s), ~ dox, contrast = c("dox", "plus", "minus"))
cat(sprintf("  mafb_dox_all: %s\n", paste(s$embryo, s$n, sep = ":", collapse = ", ")))

de_summ <- imap_dfr(res, ~ tibble(contrast = .y, genes = nrow(.x), padj_lt_0.05 = sum(.x$padj < 0.05, na.rm = TRUE),
                                  padj_lt_0.1 = sum(.x$padj < 0.1, na.rm = TRUE), pvalue_lt_0.001 = sum(.x$pvalue < 0.001, na.rm = TRUE),
                                  transgene_tf_lfc = .x$lfc[match(c(MAFB = "MAFB", EMX2 = "EMX2")[[ifelse(grepl("EMX2", .y), "EMX2", "MAFB")]], .x$gene)]))
cat("\nDE summary (min plasmid+ cells = ", MIN_POS, "):\n", sep = ""); print(as.data.frame(de_summ), row.names = FALSE)
write_csv(de_summ %>% mutate(transgene_tf_lfc = round(transgene_tf_lfc, 3)), file.path(OE_TABLES, paste0("oe_de_summary", SFX, ".csv")))
iwalk(res, ~ write_csv(.x, file.path(OE_OUT, sprintf("oe_de_%s%s.csv.gz", .y, SFX))))
top <- imap_dfr(res, ~ .x %>% filter(!is.na(padj)) %>% arrange(pvalue) %>% head(25) %>% mutate(contrast = .y)) %>%
  select(contrast, gene, baseMean, lfc, pvalue, padj) %>% mutate(across(c(baseMean, lfc), ~ round(.x, 3)), across(c(pvalue, padj), ~ signif(.x, 3)))
write_csv(top, file.path(OE_TABLES, paste0("oe_de_top25", SFX, ".csv")))

## 3. song-axis enrichment --------------------------------------------------------------------------------------------
sig <- read_csv(file.path(OE_TABLES, "reference_signature_genes.csv"), show_col_types = FALSE) %>% filter(has_ortholog, !excluded_tf)
sets <- c(sig %>% filter(specific) %>% split(paste0(.$axis, "__", .$direction, "_STRICT")) %>% map("gene"),
          sig %>% split(paste0(.$axis, "__", .$direction, "_ALL")) %>% map("gene"))
ck2zf <- oma_ck_to_zf()
gsea <- imap_dfr(res, function(r, nm) {
  st <- r %>% filter(!is.na(stat), gene %in% names(ck2zf)) %>% mutate(zf = ck2zf[gene]) %>% distinct(zf, .keep_all = TRUE)
  fg <- fgseaMultilevel(sets, setNames(st$stat, st$zf), minSize = 5, maxSize = 2000, eps = 0)
  as_tibble(fg) %>% transmute(contrast = nm, set = pathway, size, NES, pval, padj, leading_edge = map_chr(leadingEdge, ~ paste(head(.x, 15), collapse = ";")))
}) %>% separate(set, c("axis", "set"), sep = "__", remove = FALSE) %>% group_by(contrast) %>% mutate(padj = p.adjust(pval, "BH")) %>% ungroup()
cat("\nfgsea of adult song-axis signatures on each OE ranking (STRICT sets; RA-like = SONG NES > 0, SISTER NES < 0):\n")
print(as.data.frame(gsea %>% filter(grepl("STRICT", set)) %>% select(contrast, axis, set, size, NES, pval, padj) %>%
                      mutate(across(c(NES), ~ round(.x, 2)), across(c(pval, padj), ~ signif(.x, 2)))), row.names = FALSE)
write_csv(gsea %>% mutate(NES = round(NES, 3), across(c(pval, padj), ~ signif(.x, 3))), file.path(OE_TABLES, paste0("oe_fgsea_song_axes", SFX, ".csv")))
cat("\ndone\n")
