# scripts/cluster_curation.R
#
# Curation of subclusters for the reference snRNA dataset: score the things that are currently
# eyeballed (mitochondrial load, doublet load, marker co-expression, absence of unique markers),
# propose exclusions with reasons, and record the decisions in a form that survives reclustering.
#
# This dataset is the reference other datasets take their names from, so two properties matter more
# here than they would elsewhere:
#
#   1. Decisions are recorded against **cell barcodes, not cluster ids**. Reclustering renumbers
#      clusters, so a ledger keyed on cluster id silently means something else on the next run --
#      the same failure mode as souporcell's per-library `assignment`. Cluster ids appear in the
#      ledger only as provenance.
#   2. Cells are **marked, not deleted**. `curation_status` / `exclude_reason` travel with every
#      cell and the filtering happens at use time, so downstream projects can reconcile counts and
#      see what was removed and why, rather than inheriting a mystery.
#
# Companion to scripts/cluster_naming.R, which turns the surviving clusters into names.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
})

## ct_markers and friends. Sourced with the same idiom as the notebooks so there is one copy.
if (!exists("ct_markers")) {
  source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "gene_lists.R"))
}

# ---------------------------------------------------------------------------
# Class panels
# ---------------------------------------------------------------------------

#' Cell-class marker panels, derived from `ct_markers` in gene_lists.R.
#'
#' Three deviations from `ct_markers`, each deliberate. Change them here rather than in the caller.
#'
#'   - `ct_markers$Mito` is TOP2A/NUF2/MCM5/PCNA/TYMS -- **mitotic**, not mitochondrial. It is
#'     renamed `Cycling`, because this file also scores mitochondrial content (`perc_mito`) as an
#'     exclusion criterion and the two must not be confusable.
#'   - `ct_markers$Endo` is LUM/RGS5/FLI1, which is three cell types: FLI1 endothelium, RGS5
#'     pericytes, LUM fibroblast/VLMC. Split, because a reference dataset should not export one
#'     label for three populations. Set `split_vascular = FALSE` to keep the original lumping.
#'   - SYT1 is pan-neuronal rather than glutamatergic. It stays in the Glut panel (it is a true
#'     positive for those cells) but is listed in `non_specific` and excluded from specificity
#'     scoring, so a cluster does not read as Glut on SYT1 alone.
#'
#' @return list(panels, non_specific, exclusive, states)
class_panels <- function(split_vascular = TRUE) {
  panels <- list(
    Glut  = ct_markers$Glut,
    GABA  = ct_markers$GABA,
    Astro = ct_markers$Astro,
    Oligo = ct_markers$Oligo,
    OPC   = ct_markers$OPC,
    Micro = ct_markers$Micro,
    Epen  = ct_markers$Epen
  )
  if (split_vascular) {
    panels <- c(panels, list(Endo = "FLI1", Peri = "RGS5", VLMC = "LUM"))
  } else {
    panels <- c(panels, list(Endo = ct_markers$Endo))
  }

  ## States, not classes: a cluster can legitimately be Glut *and* cycling, or GABA *and* a
  ## progenitor. Scored and reported, but never counted as a co-expression violation.
  states <- list(
    Cycling = ct_markers$Mito, # see note above -- mitotic, not mitochondrial
    NP      = ct_markers$NP,
    IP      = ct_markers$IP,
    GABA_NP = ct_markers$GABA_NP
  )

  ## Sets whose members should not both be present in one cluster. Within a set, two present
  ## classes means a doublet cluster or a merge artefact. Across sets nothing is implied.
  exclusive <- list(
    neurotransmitter = c("Glut", "GABA"),
    lineage = c("Glut", "GABA", "Astro", "Oligo", "OPC", "Micro", "Epen",
                if (split_vascular) c("Endo", "Peri", "VLMC") else "Endo")
  )

  list(
    panels = panels,
    non_specific = c("SYT1"),
    exclusive = exclusive,
    states = states
  )
}

# ---------------------------------------------------------------------------
# Per-cluster QC
# ---------------------------------------------------------------------------

#' Per-cluster pseudobulk DE (one-vs-rest), plus single-cell detection fractions.
#'
#' Pseudobulk, not per-cell Wilcoxon: cells are pseudoreplicates of their individual, not
#' independent observations, so a per-cell test (this file used presto::wilcoxauc until this
#' change) treats a thousand cells from one bird as a thousand-fold more evidence than they
#' actually are and returns absurdly small p-values for anything with a real cross-individual
#' difference. Counts are summed to one profile per (cluster, individual) via
#' Seurat::AggregateExpression(), and DESeq2 -- via Seurat's own FindAllMarkers(test.use =
#' "DESeq2") wrapper, not a hand-rolled call -- tests across individuals as the actual
#' replicates.
#'
#' Three non-obvious things had to be gotten right, in order of how much they cost to find:
#'
#'  1. `slot = "counts"`, not the default `"data"`. DESeq2 fits its own negative-binomial model
#'     and needs raw counts; Seurat's DESeq2 test constructs a DESeqDataSet directly from
#'     whatever slot/layer FindMarkers hands it, so testing on log-normalized data silently
#'     produces meaningless results instead of an error. (Also note the Seurat argument is
#'     spelled `slot=`, even though `GetAssayData()` itself now requires `layer=` -- an
#'     inconsistency internal to this Seurat version, not a typo here.)
#'  2. `return.thresh = 1`. FindAllMarkers defaults to reporting only genes with p_val <= 0.01,
#'     silently dropping everything else from its return value. naming_candidates() and
#'     panel_presence() below need the *whole* per-cluster gene table -- a class-marker gene
#'     that isn't differentially expressed in some cluster must still appear as a not-present
#'     row, not vanish and read as "never tested."
#'  3. Pseudobulk pct.1/pct.2 (computed by FindAllMarkers on the aggregated object) means
#'     "fraction of pseudobulk pools with nonzero summed count", not "fraction of single cells
#'     expressing" -- coarse (moves in units of 1/n_individuals) and close to 1.0 for almost any
#'     gene once dozens of cells are summed per pool. `pct_in`/`pct_out` are therefore computed
#'     here directly from the original single-cell data, independent of the DESeq2 call, which
#'     is also why they survive unchanged from the presto-based version.
#'
#' A degenerate case worth knowing about rather than debugging blind if it recurs: synthetic
#' data with *zero* real between-replicate variance (a clean on/off spike and nothing else) can
#' trigger quasi-complete separation in DESeq2's negative-binomial fit, pinning the gene's
#' dispersion at a ceiling value and returning a large, wrong p-value for an enormous true
#' effect. Real biological replicates always carry some genuine variance, which is exactly what
#' pseudobulk testing exists to account for, so this should not arise on real data -- it is a
#' synthetic-test-only trap, not a property of the method.
#'
#' @param obj Seurat object (a subset: the comparison is always against sibling clusters)
#' @param cluster_col metadata column holding the clustering to curate
#' @param replicate_col metadata column holding the biological replicate (bird) -- pseudobulk DE
#'   has no meaning without one, so this is a hard requirement, not an optional grouping
#' @param assay assay to test in; defaults to the object's DefaultAssay. Must have a `counts`
#'   layer of raw (not normalized) counts.
#' @param min_reps clusters with fewer than this many distinct replicates error out rather than
#'   silently running an underpowered or undefined test
cluster_de <- function(obj, cluster_col, replicate_col = "individual", assay = NULL, min_reps = 2) {
  stopifnot(cluster_col %in% colnames(obj@meta.data))
  stopifnot(
    "replicate_col not found -- pseudobulk DE requires a biological-replicate column (e.g. individual)" =
      replicate_col %in% colnames(obj@meta.data)
  )
  assay <- assay %||% Seurat::DefaultAssay(obj)
  groups <- as.character(obj@meta.data[[cluster_col]])

  rep_tab <- table(groups, as.character(obj@meta.data[[replicate_col]]))
  n_reps <- rowSums(rep_tab > 0)
  too_few <- names(n_reps)[n_reps < min_reps]
  if (length(too_few) > 0) {
    stop(
      "cluster_de: fewer than ", min_reps, " '", replicate_col, "' replicates in: ",
      paste(too_few, collapse = ", "),
      " -- pseudobulk DE needs multiple biological replicates per group to estimate dispersion.",
      call. = FALSE
    )
  }

  ## Single-cell detection fractions -- see point 3 above for why these do not come from the
  ## pseudobulk object. Sparse matrix multiplication against a cluster indicator, not a
  ## split()/sapply() loop, because this runs over the full gene set.
  mat <- SeuratObject::LayerData(obj, assay = assay, layer = "data")
  detect <- mat > 0
  cluster_f <- factor(groups)
  ind <- Matrix::sparseMatrix(
    i = seq_along(cluster_f), j = as.integer(cluster_f), x = 1,
    dims = c(length(cluster_f), nlevels(cluster_f)), dimnames = list(NULL, levels(cluster_f))
  )
  n_in <- as.numeric(table(cluster_f))
  n_total <- length(cluster_f)
  sum_in <- as.matrix(detect %*% ind) # genes x clusters
  sum_total <- Matrix::rowSums(detect) # genes
  pct_in <- sweep(sum_in, 2, n_in, "/") * 100
  pct_out <- sweep(-sum_in + sum_total, 2, n_total - n_in, "/") * 100
  detection <- as_tibble(pct_in, rownames = "feature") %>%
    tidyr::pivot_longer(-feature, names_to = "cluster", values_to = "pct_in") %>%
    left_join(
      as_tibble(pct_out, rownames = "feature") %>%
        tidyr::pivot_longer(-feature, names_to = "cluster", values_to = "pct_out"),
      by = c("feature", "cluster")
    )

  ## Pseudobulk: sum raw counts within each (cluster, replicate) cell. AggregateExpression's
  ## default is already to sum from the "counts" layer (verified directly against the raw matrix);
  ## passing layer = "counts" explicitly collides with that default ("formal argument 'layer'
  ## matched by multiple actual arguments") rather than overriding it, since AggregateExpression
  ## has no `layer` formal of its own and forwards it to PseudobulkExpression() twice.
  agg <- Seurat::AggregateExpression(
    obj,
    group.by = c(cluster_col, replicate_col), assays = assay,
    return.seurat = TRUE, verbose = FALSE
  )

  ## DESeq2 refuses non-integer counts ("some values in assay are not integers"), and this is not
  ## a floating-point artefact of the summing here -- checked directly against a real dataset
  ## (Zaremba et al.) and found genuine fractional values already in its raw "counts" layer (e.g.
  ## exact .5s), the signature of upstream ambient-RNA correction or multi-mapped-read splitting.
  ## Round at the pseudobulk level, immediately before the DESeq2 call, rather than rounding the
  ## per-cell input -- standard practice for count-like data that is not strictly integer (the
  ## same accommodation DESeq2's own docs describe for RSEM/Salmon-style expected counts), and
  ## rounding after summing loses less information than rounding thousands of small per-cell
  ## fractions first.
  counts_layer <- SeuratObject::LayerData(agg, assay = assay, layer = "counts")
  SeuratObject::LayerData(agg, assay = assay, layer = "counts") <- round(counts_layer)

  fm <- suppressWarnings(Seurat::FindAllMarkers(
    agg,
    test.use = "DESeq2", assay = assay, slot = "counts", group.by = cluster_col,
    logfc.threshold = 0, min.pct = 0, return.thresh = 1, verbose = FALSE
  ))

  ## Seurat sanitizes identity-class labels containing underscores by replacing them with dashes
  ## ("Names of identity class contain underscores ('_'), replacing with dashes ('-')") -- silent
  ## except for that one-time message, and this project's cluster names commonly contain
  ## underscores (Ex_DACH2_SLIT2, cluster_glut_nido_pca_20_0.7, ...). Left unfixed, `fm$cluster`
  ## comes back dash-ified while `detection` (built from the untouched metadata column further
  ## up) keeps the real names, so the join below silently matches nothing at all -- checked
  ## directly against a real run: every single row came back with NA pct_in/pct_out. Map back to
  ## the original names via the exact substitution Seurat documents, then assert the result is a
  ## complete, unique correspondence rather than trust that assumption blindly.
  original_names <- unique(groups)
  sanitized_lookup <- setNames(original_names, gsub("_", "-", original_names))
  stopifnot(
    "cluster name sanitization is not a 1:1 mapping for this cluster_col -- Seurat's underscore/dash substitution collided two distinct names" =
      length(sanitized_lookup) == length(original_names)
  )

  fm <- fm %>%
    as_tibble() %>%
    dplyr::rename(feature = gene, logFC = avg_log2FC, padj = p_val_adj) %>%
    dplyr::mutate(cluster = unname(sanitized_lookup[as.character(cluster)])) %>%
    dplyr::select(cluster, feature, logFC, padj)

  stopifnot(
    "cluster names did not map back cleanly after undoing Seurat's underscore/dash sanitization" =
      !anyNA(fm$cluster)
  )

  fm %>% dplyr::left_join(detection, by = c("cluster", "feature"))
}

#' Per-cluster QC metrics: the quantities currently assessed by eye.
#'
#' Mitochondrial load is scored **relative to the subset**, not against a global cut. A fixed
#' threshold means different things for neurons and glia, and an upstream perc_mito filter has
#' already run, so what is left to detect is enrichment within this subset. `mito_mad` is
#' (cluster median - subset median) / MAD; MAD is used rather than SD because a bad cluster would
#' inflate SD and hide itself.
#'
#' @param panels output of class_panels()
#' @param pct_present fraction of a cluster's cells that must express a panel's genes (median over
#'   the panel) for that class to count as present
cluster_qc <- function(obj, cluster_col, de = NULL, panels = class_panels(),
                       mito_col = "perc_mito", pct_present = 0.25, assay = NULL) {
  de <- de %||% cluster_de(obj, cluster_col, assay = assay)
  md <- obj@meta.data %>%
    as_tibble(rownames = "cell") %>%
    mutate(cluster = as.character(.data[[cluster_col]]))

  ## presto reports pct_in on a 0-100 scale.
  present <- panel_presence(de, panels, pct_present = pct_present * 100)

  base <- md %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      median_mito = if (mito_col %in% names(md)) median(.data[[mito_col]], na.rm = TRUE) else NA_real_,
      median_nFeature = if ("nFeature_RNA" %in% names(md)) median(nFeature_RNA, na.rm = TRUE) else NA_real_,
      ## Two independent doublet calls: DoubletFinder (expression-based) and souporcell
      ## (genotype-based). By the time any subset object exists these are **expected to be zero** --
      ## snrna_seurat_cellbender_preprocess.qmd keeps only the intersection of the two callers'
      ## singlets before subclustering. They are reported as a *verification* that the upstream
      ## filter ran, not as a curation criterion: a non-zero value here means the object did not go
      ## through that step, and the flagging thresholds below deliberately do not act on them.
      ## Residual doublets that survive both callers are what the marker co-expression rule catches.
      frac_df_doublet = if ("DF_class" %in% names(md)) mean(DF_class == "Doublet", na.rm = TRUE) else NA_real_,
      frac_soc_doublet = if ("doublet_status" %in% names(md)) mean(doublet_status == "doublet", na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  mito_med <- median(base$median_mito, na.rm = TRUE)
  mito_mad <- stats::mad(base$median_mito, na.rm = TRUE)

  base %>%
    mutate(
      ## MAD of 0 (every cluster identical) would divide to Inf; report NA rather than flag everything.
      mito_mad = if (is.na(mito_mad) || mito_mad == 0) NA_real_ else (median_mito - mito_med) / mito_mad
    ) %>%
    left_join(present, by = "cluster") %>%
    arrange(dplyr::desc(n_cells))
}

#' Which class panels are "present" in each cluster, and whether that violates exclusivity.
panel_presence <- function(de, panels = class_panels(), pct_present = 25) {
  specific <- map(panels$panels, ~ setdiff(.x, panels$non_specific))
  specific <- specific[map_int(specific, length) > 0]

  score_one <- function(genes) {
    de %>%
      dplyr::filter(feature %in% genes) %>%
      group_by(cluster) %>%
      summarise(pct = median(pct_in, na.rm = TRUE), .groups = "drop")
  }

  scores <- imap(specific, ~ score_one(.x) %>% mutate(class = .y)) %>% bind_rows()
  state_scores <- imap(panels$states, ~ score_one(.x) %>% mutate(class = .y)) %>% bind_rows()

  wide <- scores %>%
    mutate(present = pct >= pct_present) %>%
    group_by(cluster) %>%
    summarise(
      classes_present = paste(sort(class[present]), collapse = "+"),
      n_classes_present = sum(present),
      top_class = class[which.max(pct)],
      top_class_pct = round(max(pct), 1),
      ## Margin over the runner-up: a class call with no margin is a coin flip, and this is the
      ## number to look at before believing `top_class`.
      class_margin = round(max(pct) - sort(pct, decreasing = TRUE)[2], 1),
      .groups = "drop"
    )

  violations <- scores %>%
    dplyr::filter(pct >= pct_present) %>%
    group_by(cluster) %>%
    summarise(
      coexpr_violation = any(map_lgl(panels$exclusive, ~ sum(class %in% .x) > 1)),
      .groups = "drop"
    )

  states_wide <- state_scores %>%
    dplyr::filter(pct >= pct_present) %>%
    group_by(cluster) %>%
    summarise(states_present = paste(sort(class), collapse = "+"), .groups = "drop")

  wide %>%
    left_join(violations, by = "cluster") %>%
    left_join(states_wide, by = "cluster") %>%
    mutate(
      coexpr_violation = coalesce(coexpr_violation, FALSE),
      states_present = coalesce(states_present, "")
    )
}

# ---------------------------------------------------------------------------
# Flagging
# ---------------------------------------------------------------------------

#' Default flagging thresholds. Passed explicitly so a run records what it used.
#'
#' `doublet_frac` defaults to NA, i.e. **disabled**, because the preprocessing notebook already
#' subsets to cells both DoubletFinder and souporcell call singlets before any subclustering
#' happens. Every subset object therefore has a doublet fraction of exactly zero in every cluster,
#' and a rule on it would be a check that can never fire -- worse than no check, because it reads
#' like coverage. Set it to a number only when curating an object that has *not* been through that
#' filter. Doublets that survive both callers show up as marker co-expression, which is a separate
#' rule and is not disabled.
curation_thresholds <- function(mito_mad = 3, mito_abs = NA_real_, doublet_frac = NA_real_,
                                min_cells = 20, require_unique_markers = TRUE) {
  list(mito_mad = mito_mad, mito_abs = mito_abs, doublet_frac = doublet_frac,
       min_cells = min_cells, require_unique_markers = require_unique_markers)
}

#' Flag clusters for exclusion, with a reason per flag.
#'
#' `flag_reason` is a "+"-joined string so a cluster tripping several criteria says so; a cluster
#' flagged only by one weak criterion is the one to look at by hand. Nothing is dropped here --
#' this proposes, the ledger decides.
#'
#' @param naming optional output of select_naming_genes(); used for the no-unique-markers rule,
#'   which says *merge or drop*, not silently keep.
flag_clusters <- function(qc, thresholds = curation_thresholds(), naming = NULL) {
  out <- qc %>%
    mutate(
      flag_mito = !is.na(mito_mad) & mito_mad > thresholds$mito_mad |
        (!is.na(thresholds$mito_abs) & !is.na(median_mito) & median_mito > thresholds$mito_abs),
      ## NA threshold = rule disabled (the default; see curation_thresholds()).
      flag_doublet = !is.na(thresholds$doublet_frac) &
        (coalesce(frac_df_doublet, 0) > thresholds$doublet_frac |
           coalesce(frac_soc_doublet, 0) > thresholds$doublet_frac),
      flag_coexpr = coexpr_violation,
      flag_small = n_cells < thresholds$min_cells
    )

  if (!is.null(naming) && thresholds$require_unique_markers) {
    uniq <- naming %>% dplyr::select(cluster, n_naming_genes)
    out <- out %>%
      left_join(uniq, by = "cluster") %>%
      mutate(flag_no_markers = coalesce(n_naming_genes, 0L) == 0L)
  } else {
    out <- out %>% mutate(flag_no_markers = FALSE)
  }

  out %>%
    mutate(
      flag_reason = pmap_chr(
        list(flag_mito, flag_doublet, flag_coexpr, flag_small, flag_no_markers),
        function(m, d, c, s, k) {
          paste(c("mito_high", "doublet_high", "marker_coexpression", "too_few_cells",
                  "no_unique_markers")[c(m, d, c, s, k)], collapse = "+")
        }
      ),
      flagged = flag_reason != ""
    )
}

# ---------------------------------------------------------------------------
# The ledger
# ---------------------------------------------------------------------------

LEDGER_COLS <- c("cell", "subset", "round", "cluster_at_decision", "status", "reason", "decided_on")

#' Expand flagged clusters into a per-cell exclusion proposal.
#'
#' The expansion to barcodes is the point: it is what makes the decision survive the reclustering
#' that follows it.
propose_exclusions <- function(obj, cluster_col, flagged, subset_name, round,
                               decided_on = Sys.Date()) {
  drop <- flagged %>% dplyr::filter(flagged)
  md <- obj@meta.data %>%
    as_tibble(rownames = "cell") %>%
    mutate(cluster = as.character(.data[[cluster_col]]))

  md %>%
    inner_join(drop %>% dplyr::select(cluster, reason = flag_reason), by = "cluster") %>%
    transmute(
      cell, subset = subset_name, round = round, cluster_at_decision = cluster,
      status = "exclude", reason, decided_on = as.character(decided_on)
    )
}

read_ledger <- function(path) {
  if (!file.exists(path)) {
    return(tibble(!!!set_names(rep(list(character()), length(LEDGER_COLS)), LEDGER_COLS)))
  }
  readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()))
}

#' Append decisions, keeping one row per (cell, subset) -- the latest round wins.
write_ledger <- function(ledger, new_rows, path) {
  combined <- bind_rows(ledger, mutate(new_rows, across(everything(), as.character))) %>%
    group_by(cell, subset) %>%
    slice_max(as.integer(round), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(subset, as.integer(round), cell)
  readr::write_csv(combined, path)
  combined
}

#' Stamp `curation_status` / `exclude_reason` onto an object from the ledger.
#'
#' Marks rather than subsets. Cells absent from the ledger are "keep" -- a ledger records
#' exceptions, so a missing row is not a missing decision.
apply_ledger <- function(obj, ledger, subset_name = NULL) {
  led <- if (is.null(subset_name)) ledger else dplyr::filter(ledger, subset == subset_name)
  status <- setNames(led$status, led$cell)
  reason <- setNames(led$reason, led$cell)
  cells <- colnames(obj)
  obj$curation_status <- unname(coalesce(status[cells], "keep"))
  obj$exclude_reason <- unname(coalesce(reason[cells], NA_character_))
  obj
}

#' Cells that survive curation.
kept_cells <- function(obj) colnames(obj)[obj$curation_status == "keep"]

# ---------------------------------------------------------------------------
# Did the curation delete biology?
# ---------------------------------------------------------------------------

#' Composition of the excluded cells against the design.
#'
#' The check that matters for this dataset. `lman`/`nr` are batch 2 -- different birds, different
#' ambient and doublet rates -- so a quality-driven rule can preferentially delete batch-2 cells,
#' and the result would read as "LMAN lacks cell type X". If exclusions concentrate in one batch,
#' bird or region, that is a finding to report, not a filter to apply.
#'
#' @param vars metadata columns to test composition against
excluded_composition <- function(obj, vars = c("batch", "individual", "position", "sample_name")) {
  vars <- intersect(vars, colnames(obj@meta.data))
  md <- obj@meta.data %>% as_tibble(rownames = "cell")
  overall <- mean(md$curation_status == "exclude")

  map(vars, function(v) {
    md %>%
      group_by(level = as.character(.data[[v]])) %>%
      summarise(
        n_cells = n(),
        n_excluded = sum(curation_status == "exclude"),
        frac_excluded = n_excluded / n_cells,
        .groups = "drop"
      ) %>%
      mutate(
        variable = v,
        overall_frac = overall,
        ## Ratio to the overall rate: 1 means this level is losing cells at the average rate.
        enrichment = ifelse(overall > 0, frac_excluded / overall, NA_real_)
      ) %>%
      dplyr::select(variable, level, n_cells, n_excluded, frac_excluded, overall_frac, enrichment)
  }) %>%
    bind_rows() %>%
    arrange(dplyr::desc(enrichment))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
