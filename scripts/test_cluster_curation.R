## Synthetic test of scripts/cluster_curation.R + scripts/cluster_naming.R.
## Structure is planted, so the right answer is known in advance.
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(tibble); library(purrr)})
setwd("/ssd/brad/rstudio/multiome/song-system-grn")
source("scripts/cluster_curation.R")
source("scripts/cluster_naming.R")

set.seed(42)
n_per <- 120
clusters <- c("A", "B", "C", "D")
genes <- c("SLC17A6", "SYT1", "GAD1", "GAD2", "DLX1", "SLC1A3", "GFAP", "MBP", "PLP1",
           "PDGFRA", "CSF1R", "FLI1", "RGS5", "LUM", "SPEF2", "FOXJ1", "TOP2A", "PCNA",
           "NUF2", "MCM5", "TYMS", "SOX2", "PAX6", "NES", "CCND2", "CDCA7", "VAX1",
           "NKX2-1", "OLIG2", "ISL1", "EOMES",
           paste0("GENE", 1:40), "SHARED1", "UNIQA", "UNIQB", "PAIR1", "PAIR2",
           ## Pure background padding, tested by nothing below: DESeq2's dispersion trend fit
           ## (fitType = "local", a locfit curve over mean expression) needs enough genes to be
           ## numerically stable. ~80 genes is not enough -- it either crashes ("newsplit: out of
           ## vertex space") or fits a degenerate trend; real single-cell data has thousands of
           ## genes and never hits this, so this padding exists only to make a small synthetic
           ## fixture behave like real data for the test that follows.
           paste0("BGPAD", 1:1500))
cells <- paste0("cell", seq_len(n_per * length(clusters)))
lab <- rep(clusters, each = n_per)
individual <- rep(rep(paste0("bird", 1:3), each = 40), 4)

## Per-(gene, bird) multiplicative noise (log-normal, ~30% CV): real biological replicate
## variance. Without it, a clean deterministic on/off spike causes quasi-complete separation in
## DESeq2's negative-binomial fit, pinning dispersion at a ceiling and returning a large p-value
## for an enormous true effect -- a synthetic-data trap documented in cluster_de(), and this
## fixture must not reproduce it by being unrealistically noise-free.
birds <- unique(individual)
bird_effect <- matrix(rlnorm(length(genes) * length(birds), 0, 0.3), nrow = length(genes),
                      dimnames = list(genes, birds))

base_rate <- matrix(0.05, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
on <- function(g, cl, lambda = 25) {
  idx <- which(lab %in% cl)
  base_rate[g, idx] <<- lambda
}
## All four are glutamatergic; D is a Glut/GABA doublet cluster (co-expression violation).
on(c("SLC17A6", "SYT1"), clusters)
on(c("GAD1", "GAD2", "DLX1"), "D")
## Naming structure: A has a unique gene; B has a unique gene; C has none of its own, only a gene
## shared with A -- so C should need a *pair* to be identified, or fail.
on("UNIQA", "A"); on("UNIQB", "B")
on("SHARED1", c("A", "C"))
on("PAIR1", c("C", "D")); on("PAIR2", c("C", "B"))   # C = PAIR1 & PAIR2, unique only in combination

rate <- base_rate * bird_effect[, individual]
counts <- matrix(rpois(length(rate), rate), nrow = nrow(rate), dimnames = dimnames(rate))

obj <- CreateSeuratObject(counts = counts)
obj$seurat_clusters_test <- lab
obj$perc_mito <- ifelse(lab == "B", rnorm(length(lab), 6, .4), rnorm(length(lab), 1, .3))  # B = high mito
obj$DF_class <- ifelse(lab == "D" & runif(length(lab)) < .6, "Doublet", "Singlet")
obj$doublet_status <- "singlet"
obj$individual <- individual
obj$position <- ifelse(lab == "A", "ra", sample(c("ra", "hvc", "nc"), length(lab), TRUE))
obj$batch <- 1L
obj <- NormalizeData(obj, verbose = FALSE)

de <- cluster_de(obj, "seurat_clusters_test")
qc <- cluster_qc(obj, "seurat_clusters_test", de = de)
nm <- select_naming_genes(de)
fl <- flag_clusters(qc, naming = nm)

cat("=== QC ===\n");     print(as.data.frame(qc %>% select(cluster, n_cells, median_mito, mito_mad, frac_df_doublet, classes_present, n_classes_present, coexpr_violation)))
cat("\n=== naming genes ===\n"); print(as.data.frame(nm))
cat("\n=== flags ===\n");  print(as.data.frame(fl %>% select(cluster, flagged, flag_reason)))

cat("\n=== ASSERTIONS ===\n")
chk <- function(lab, ok) cat(sprintf("%-58s %s\n", lab, ifelse(ok, "PASS", "FAIL")))
chk("B flagged for high mito", grepl("mito_high", fl$flag_reason[fl$cluster == "B"]))
chk("D flagged for Glut/GABA co-expression", grepl("marker_coexpression", fl$flag_reason[fl$cluster == "D"]))
## The doublet rule is disabled by default, because every subset object reaching curation has
## already been filtered to both callers' singlets upstream -- a rule that can never fire reads as
## coverage it does not provide. Assert both halves: off by default, and working when asked for.
chk("doublet rule OFF by default", !grepl("doublet_high", fl$flag_reason[fl$cluster == "D"]))
fl_dbl <- flag_clusters(qc, thresholds = curation_thresholds(doublet_frac = 0.25), naming = nm)
chk("doublet rule fires when enabled", grepl("doublet_high", fl_dbl$flag_reason[fl_dbl$cluster == "D"]))
chk("D still caught by co-expression with the rule off",
    grepl("marker_coexpression", fl$flag_reason[fl$cluster == "D"]))
chk("A and C not flagged for mito", !any(grepl("mito", fl$flag_reason[fl$cluster %in% c("A","C")])))
chk("A named by its unique gene", grepl("UNIQA", nm$naming_genes[nm$cluster == "A"]))
chk("B named by its unique gene", grepl("UNIQB", nm$naming_genes[nm$cluster == "B"]))
chk("A uniquely identified by 1 gene", nm$n_naming_genes[nm$cluster == "A"] == 1 && nm$unique_id[nm$cluster == "A"])
## C's only candidate markers (SHARED1, PAIR1, PAIR2) are each shared with exactly one other
## cluster (A, D, B respectively) -- under presto's rank-based one-vs-rest test this combined
## into a unique 2-gene name; under one-vs-rest pseudobulk DESeq2 it legitimately does not.
## Mechanism, confirmed directly against the DE table: at cluster C, PAIR1 shows logFC=1.69,
## pct_in=100, but padj=1 -- because PAIR1 is *also* genuinely elevated in D, so the "rest"
## pool (A+B+D) it is tested against is a real mixture of low and high values, which inflates
## the working dispersion for that one-vs-rest contrast far more for a parametric NB test than
## for a rank statistic. The same gene is highly significant (padj as low as 8e-7) at the
## clusters where it truly is private. This is structural to combining one-vs-rest testing with
## combination-naming, not specific to these genes or fixable by loosening naming_filters() --
## any "half-unique, shared with one specific sibling" marker hits it. Documented here rather
## than propped up by retuning thresholds to one synthetic case.
chk("C gets no name under one-vs-rest pseudobulk DESeq2 (structural, not a bug)",
    nm$n_naming_genes[nm$cluster == "C"] == 0 && !nm$unique_id[nm$cluster == "C"])
chk("C is flagged no_unique_markers as a result",
    grepl("no_unique_markers", fl$flag_reason[fl$cluster == "C"]))
chk("SYT1 never used as a naming gene (pan-neuronal)", !any(grepl("SYT1", nm$naming_genes), na.rm = TRUE))
chk("Cycling panel is mitotic, not perc_mito", identical(class_panels()$states$Cycling, ct_markers$Mito))
chk("vascular split into Endo/Peri/VLMC", all(c("Endo","Peri","VLMC") %in% names(class_panels()$panels)))

## Ledger round-trip: decisions must survive renumbering of clusters. Scoped to mito/co-expression
## flags specifically (require_unique_markers = FALSE) -- this section tests ledger *mechanics*
## (does an exclusion decision survive relabeling), not naming quality, and C's no_unique_markers
## flag is a real, DE-engine-dependent outcome tested separately above; coupling the two here would
## make ledger-mechanics assertions fragile to changes in naming power that have nothing to do with
## the ledger.
fl_ledger <- flag_clusters(qc, thresholds = curation_thresholds(require_unique_markers = FALSE))
led_path <- tempfile(fileext = ".csv")
prop <- propose_exclusions(obj, "seurat_clusters_test", fl_ledger, "test_subset", round = 1)
led <- write_ledger(read_ledger(led_path), prop, led_path)
obj2 <- apply_ledger(obj, led, "test_subset")
excluded_cells <- colnames(obj2)[obj2$curation_status == "exclude"]

## Now renumber every cluster and re-apply: the same cells must be excluded.
obj3 <- obj
obj3$seurat_clusters_test <- paste0("X", match(obj$seurat_clusters_test, clusters) * 7)
obj3 <- apply_ledger(obj3, read_ledger(led_path), "test_subset")
chk("ledger survives cluster renumbering (same cells)",
    setequal(excluded_cells, colnames(obj3)[obj3$curation_status == "exclude"]))
chk("excluded cells are B and D only",
    setequal(unique(lab[match(excluded_cells, cells)]), c("B", "D")))
chk("kept_cells() excludes them", length(kept_cells(obj2)) == length(cells) - length(excluded_cells))

comp <- excluded_composition(obj2)
chk("excluded_composition reports per-variable enrichment", nrow(comp) > 0 && "enrichment" %in% names(comp))
cat("\n=== excluded composition (top) ===\n"); print(as.data.frame(head(comp, 6)))

## --- Zaremba-style alias: additive, non-destructive, overrides win over generated names ---
zl <- zaremba_labels(obj, "seurat_clusters_test", de = de,
                     overrides = c(A = "Ex_XENIUM_OVERRIDE"), override_source = "xenium (synthetic)")
cat("\n=== Zaremba-style aliases ===\n"); print(as.data.frame(zl))

chk("zaremba_style prefixes Ex_ and joins with underscore",
    zaremba_style("UNIQA/UNIQB") == "Ex_UNIQA_UNIQB")
chk("zaremba_style caps at max_genes (2 of 3)",
    zaremba_style("G1/G2/G3", max_genes = 2) == "Ex_G1_G2")
chk("zaremba_style NA-safe for unnamed clusters", is.na(zaremba_style(NA_character_)))
chk("override wins over generated name for A", zl$celltype_ex_style[zl$cluster == "A"] == "Ex_XENIUM_OVERRIDE")
chk("override source recorded for A", zl$source[zl$cluster == "A"] == "xenium (synthetic)")
chk("B gets a generated name, not an override",
    zl$source[zl$cluster == "B"] == "generated from own markers" &&
      zl$celltype_ex_style[zl$cluster == "B"] == "Ex_UNIQB")
chk("un-overridden, unnamed cluster reports its own reason",
    all(is.na(zl$celltype_ex_style) == (zl$source == "no marker cleared naming_filters()")))
chk("original cluster labels A-D untouched by aliasing",
    setequal(zl$cluster, clusters))

## overrides = NULL is what a bare c() placeholder evaluates to -- a very natural way to write
## "no overrides yet" -- and previously broke coalesce() with "object 'overridden' not found"
## instead of behaving like character(0). Caught rendering glut_ex_style_naming.qmd for real.
zl_null <- tryCatch(
  zaremba_labels(obj, "seurat_clusters_test", de = de, overrides = NULL),
  error = function(e) e
)
chk("zaremba_labels(overrides = NULL) does not error", !inherits(zl_null, "error"))
chk("zaremba_labels(overrides = NULL) behaves like no overrides",
    !inherits(zl_null, "error") && identical(
      zl_null$celltype_ex_style,
      zaremba_labels(obj, "seurat_clusters_test", de = de, overrides = character())$celltype_ex_style
    ))
