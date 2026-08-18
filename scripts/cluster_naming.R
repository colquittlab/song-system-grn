# scripts/cluster_naming.R
#
# Naming clusters from their own differential expression, in the form papers actually use: take the
# DEGs across a subset, then keep the 1-3 genes that *uniquely* identify each cluster.
#
# The word doing the work is "uniquely". A gene being the top DEG of a cluster is not enough -- the
# same gene is often top in a sibling too, and a name built on it does not distinguish them. So
# selection here is a small set-cover problem: choose the fewest genes whose joint presence is true
# of this cluster and of no sibling. That is why a cluster sometimes gets a two- or three-gene name
# and its neighbour gets one; the length of the name is a measurement, not a style choice.
#
# Companion to scripts/cluster_curation.R. Both consume one presto pass (cluster_de()).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
})

#' Filters for what may serve as a naming gene. Passed explicitly so a run records its choices.
#'
#' `max_present_clusters` is the constraint that makes combination naming possible, and getting it
#' wrong quietly disables this whole file. The obvious formulation -- a low `pct_out_max`, "the gene
#' must be off in the rest of the subset" -- rejects *every* gene shared with any sibling, because a
#' gene shared with one of four clusters already sits at pct_out ~33. That leaves only genes that
#' are private to one cluster, so the set-cover step below can never combine anything and any
#' cluster without a private marker gets no name at all. (Caught exactly this way in testing: a
#' planted cluster identifiable only by PAIR1+PAIR2 came back unnamed.)
#'
#' The criterion that expresses the intent is a cap on the number of *clusters* a gene is present
#' in: a gene on in 2 of 12 clusters is a fine half of a two-gene name, while a gene on in 9 of 12
#' names nothing. `pct_out_max` stays as a loose backstop against ubiquitous genes.
#'
#' @param max_present_clusters absolute cap; NULL means ceiling(n_clusters / 3), min 2
naming_filters <- function(pct_in_min = 40, pct_out_max = 60,
                           logFC_min = 0.25, padj_max = 0.05, max_genes = 3,
                           max_present_clusters = NULL) {
  list(pct_in_min = pct_in_min, pct_out_max = pct_out_max,
       logFC_min = logFC_min, padj_max = padj_max, max_genes = max_genes,
       max_present_clusters = max_present_clusters)
}

#' Candidate naming genes per cluster, ranked.
#'
#' Rank is significance (padj) first, then the pct_in/pct_out gap. There is no AUC here -- `de`
#' comes from cluster_de()'s pseudobulk DESeq2 test (see that function), which reports a p-value
#' and a log2 fold change, not a rank statistic. `logFC_min`/`padj_max` were already presto-scale
#' thresholds carried over from that test; they are ordinary values for a DESeq2 log2FoldChange
#' and BH-adjusted p-value too, so they are kept as the defaults, but note the *meaning* of
#' "0.25" changed with the test -- it is now a log2 fold change on pseudobulk sums, not a
#' natural-log difference of single-cell means.
naming_candidates <- function(de, filters = naming_filters(), exclude_genes = character(),
                              pct_present = 25) {
  n_clusters <- dplyr::n_distinct(de$cluster)
  cap <- filters$max_present_clusters %||% max(2L, ceiling(n_clusters / 3))

  ## How many clusters is each gene on in? Computed on the same presence rule the uniqueness step
  ## uses, so candidacy and uniqueness are judged on one definition of "expressed".
  n_present <- de %>%
    group_by(feature) %>%
    summarise(n_clusters_present = sum(pct_in >= pct_present), .groups = "drop")

  de %>%
    left_join(n_present, by = "feature") %>%
    dplyr::filter(
      !feature %in% exclude_genes,
      !is.na(padj), # DESeq2 sets padj to NA for genes it excludes via independent filtering
      pct_in >= filters$pct_in_min,
      pct_out <= filters$pct_out_max,
      logFC >= filters$logFC_min,
      padj <= filters$padj_max,
      n_clusters_present <= cap
    ) %>%
    ## Rank prefers genes on in fewer clusters, then significance, then the pct_in/pct_out gap: a
    ## gene private to this cluster should beat one shared with three others even at equal padj.
    mutate(gap = pct_in - pct_out, score = -log10(pmax(padj, 1e-300)) + gap / 200) %>%
    arrange(cluster, n_clusters_present, dplyr::desc(score)) %>%
    dplyr::select(cluster, feature, logFC, pct_in, pct_out, gap, score, padj, n_clusters_present)
}

#' Logical gene x cluster matrix: is this gene "on" in this cluster?
#'
#' Deliberately a *lower* bar than the candidate filter. Candidacy asks "could this gene name this
#' cluster"; presence asks "would a reader looking at a violin plot call this gene expressed here".
#' Uniqueness has to be judged against the second, or a name looks unique only because the sibling
#' narrowly missed a significance cut.
presence_matrix <- function(de, pct_present = 25) {
  de %>%
    mutate(present = pct_in >= pct_present) %>%
    dplyr::select(feature, cluster, present) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = present, values_fill = FALSE)
}

#' Choose the smallest set of candidate genes that identifies each cluster uniquely.
#'
#' Greedy set cover: at each step take the candidate that most shrinks the set of sibling clusters
#' still consistent with the genes chosen so far, breaking ties by rank. Greedy is not guaranteed
#' optimal, but with max_genes = 3 the search space is small and the failure mode is a name one gene
#' longer than necessary, which is harmless. Stops early when the cluster is uniquely identified.
#'
#' Returns one row per cluster with the chosen genes, how many clusters remain ambiguous, and
#' whether uniqueness was reached. `n_naming_genes == 0` means nothing passed the filters -- that
#' cluster has no marker of its own and is a merge-or-drop candidate, which is what
#' flag_clusters(require_unique_markers = TRUE) acts on.
select_naming_genes <- function(de, filters = naming_filters(), pct_present = 25,
                                exclude_genes = character()) {
  cands <- naming_candidates(de, filters, exclude_genes, pct_present = pct_present)
  pres <- presence_matrix(de, pct_present)
  clusters <- sort(unique(as.character(de$cluster)))

  present_in <- function(gene) {
    row <- pres[pres$feature == gene, , drop = FALSE]
    if (nrow(row) == 0) return(character())
    clusters[unlist(row[1, clusters]) %in% TRUE]
  }

  map_dfr(clusters, function(cl) {
    cc <- cands %>% dplyr::filter(cluster == cl)
    if (nrow(cc) == 0) {
      return(tibble(cluster = cl, naming_genes = NA_character_, n_naming_genes = 0L,
                    n_ambiguous = length(clusters), unique_id = FALSE,
                    top_padj = NA_real_, ambiguous_with = paste(setdiff(clusters, cl), collapse = ",")))
    }
    chosen <- character()
    ambiguous <- clusters
    repeat {
      remaining <- setdiff(cc$feature, chosen)
      if (length(remaining) == 0 || length(chosen) >= filters$max_genes) break
      ## Score each remaining candidate by how small it would make the ambiguity set.
      sizes <- map_int(remaining, function(g) length(intersect(ambiguous, present_in(g))))
      best <- remaining[order(sizes, match(remaining, cc$feature))][1]
      new_amb <- intersect(ambiguous, present_in(best))
      ## A gene that removes no ambiguity adds nothing to the name; stop rather than pad.
      if (length(chosen) > 0 && length(new_amb) == length(ambiguous)) break
      chosen <- c(chosen, best)
      ambiguous <- new_amb
      if (length(ambiguous) <= 1) break
    }
    tibble(
      cluster = cl,
      naming_genes = paste(chosen, collapse = "/"),
      n_naming_genes = length(chosen),
      n_ambiguous = length(ambiguous),
      unique_id = length(ambiguous) <= 1,
      top_padj = cc$padj[1],
      ambiguous_with = paste(setdiff(ambiguous, cl), collapse = ",")
    )
  })
}

#' Format naming genes in the style of Zaremba et al. 2025 (Science 387:adp5182), the chicken
#' pallium atlas -- `Ex_<gene1>[_<gene2>]` for excitatory/glutamatergic types.
#'
#' Confirmed directly from their object (`/hdd/sc_datasets/zaremba_gg_dev_2024/
#' Gg_adult_snRNA_seq_srt.rds`, metadata columns `anno_level_1`/`anno_level_3`), not from the
#' published text, which is paywalled: `Ex_` prefix, genes joined by `_`, at most two genes.
#' Their finer-level names are **not** a strict extension of the coarser one -- some supertypes
#' drop the parent subclass's gene entirely (e.g. `Ex_TSHZ2_NR4A2` shares nothing with any
#' `anno_level_1` subclass name) -- so this format function only reformats *this* dataset's own
#' `select_naming_genes()` output; it makes no claim that the result denotes the same cell type
#' Zaremba's identically-formatted name does. That claim is a cross-species homology question,
#' answered by integration or by independent (e.g. spatial) evidence, not by string formatting --
#' see `zaremba_overrides` in make_labels().
#'
#' `max_genes = 2` here, not the package default of 3: Zaremba's names never exceed two, and this
#' dataset's own resolution is not directly comparable to theirs, so matching their cap is a
#' convention choice, not a claim of equivalent granularity.
zaremba_style <- function(naming_genes, max_genes = 2) {
  ifelse(
    is.na(naming_genes) | naming_genes == "",
    NA_character_,
    paste0("Ex_", map_chr(str_split(naming_genes, "/"), ~ paste(head(.x, max_genes), collapse = "_")))
  )
}

#' Ex_-style alias per existing cluster: your own confirmed correspondences where you have them,
#' a fresh name from this dataset's own top markers everywhere else.
#'
#' Deliberately additive and non-destructive: it never touches `cluster`/`celltype` (or whatever
#' column already names the existing clusters) and never reclusters anything. It reads an existing
#' grouping, computes each group's own naming genes, and returns a lookup from that group to an
#' `Ex_`-style alias plus how the alias was decided -- for a `left_join` onto the object as a new
#' column, per the point of adopting this nomenclature without disturbing the clusters that are
#' already load-bearing for the manuscript.
#'
#' @param obj Seurat object, unmodified group of interest already present as a metadata column
#' @param cluster_col existing clustering to alias -- NOT recomputed
#' @param overrides named character vector: existing cluster label -> Zaremba name, for clusters
#'   where you already have a confident correspondence (Xenium spatial mapping, prior integration,
#'   etc.) that this function cannot derive from expression alone. Source of the correspondence is
#'   recorded as free text so the ledger says where a name came from, e.g.
#'   `c("Glut-HVC-1" = "Ex_DACH2_SLIT2")`, `override_source = "xenium 2026-08 mapping"`.
#' @param de optional precomputed cluster_de() output, to avoid rerunning presto if it already ran
#'   for curation on the same grouping.
zaremba_labels <- function(obj, cluster_col, overrides = character(), override_source = NA_character_,
                          filters = naming_filters(max_genes = 2), de = NULL, assay = NULL) {
  de <- de %||% cluster_de(obj, cluster_col, assay = assay)
  nm <- select_naming_genes(de, filters = filters)

  ## `overrides[cluster]` on a NULL indexes to NULL, and mutate() silently *drops* a column
  ## assigned NULL rather than creating an all-NA one -- so `overrides = NULL` (not this function's
  ## own default, but exactly what a bare `c()` placeholder for "no overrides yet" evaluates to,
  ## the natural way to write one before you have any) breaks `coalesce(overridden, ...)` below with
  ## "object 'overridden' not found" instead of behaving like no overrides at all. Guard here rather
  ## than trust every caller to write `character()` instead of `c()`.
  if (is.null(overrides)) overrides <- character()

  nm %>%
    mutate(
      generated = zaremba_style(naming_genes, max_genes = filters$max_genes),
      overridden = overrides[cluster],
      celltype_ex_style = coalesce(overridden, generated),
      source = case_when(
        !is.na(overridden) ~ coalesce(override_source, "manual override"),
        !is.na(generated) ~ "generated from own markers",
        TRUE ~ "no marker cleared naming_filters()"
      )
    ) %>%
    dplyr::select(cluster, celltype_ex_style, source, naming_genes, n_naming_genes, unique_id)
}

#' Assemble labels from the three fields, each of which came from its own evidence.
#'
#' `<class>-<anatomy>-<genes>`, dropping any field that was not supported. Fields are never
#' invented to fill the pattern: a cluster whose class call had no margin gets "Unresolved", and a
#' cluster with no anatomical enrichment simply has no anatomy field.
make_labels <- function(naming, qc = NULL, anatomy = NULL, class_margin_min = 10,
                        prefix = NULL) {
  out <- naming %>% dplyr::select(cluster, naming_genes, n_naming_genes, unique_id)

  if (!is.null(qc)) {
    out <- out %>%
      left_join(qc %>% dplyr::select(cluster, top_class, top_class_pct, class_margin), by = "cluster") %>%
      mutate(class = if_else(!is.na(class_margin) & class_margin >= class_margin_min,
                             top_class, "Unresolved"))
  } else {
    out <- out %>% mutate(class = NA_character_)
  }

  if (!is.null(anatomy)) {
    out <- out %>% left_join(anatomy %>% dplyr::select(cluster, region), by = "cluster")
  } else {
    out <- out %>% mutate(region = NA_character_)
  }

  out %>%
    mutate(
      label = pmap_chr(list(class, region, naming_genes), function(cl, rg, gn) {
        paste(c(prefix, cl, rg, gn)[!is.na(c(prefix, cl, rg, gn)) &
                                      nzchar(c(prefix %||% "", cl %||% "", rg %||% "", gn %||% ""))],
              collapse = "-")
      })
    )
}

# ---------------------------------------------------------------------------
# Anatomy: an enrichment call, not a marker call
# ---------------------------------------------------------------------------

#' Which region (if any) a cluster is enriched for, using birds as the replicate unit.
#'
#' Cells are not independent, so a cell-level test calls everything significant. The unit here is
#' the individual: a region label is earned when the cluster's cells come disproportionately from
#' that region *in most of the birds that contribute to it*.
#'
#' **Read the `confounded` column before using the result.** Region and batch are partly confounded
#' in this dataset -- `lman`/`nr` are batch 2 and share no birds with `ra`/`arco`/`hvc`/`nc` -- so an
#' enrichment that lives entirely in one batch cannot be distinguished from a batch effect. Those
#' rows get region = NA and confounded = TRUE rather than a label that would overstate the evidence.
region_enrichment <- function(obj, cluster_col, region_col = "position",
                              individual_col = "individual", batch_col = "batch",
                              frac_min = 0.5, min_birds = 2) {
  md <- obj@meta.data %>%
    as_tibble(rownames = "cell") %>%
    mutate(cluster = as.character(.data[[cluster_col]]))
  if (!region_col %in% names(md)) return(tibble(cluster = character(), region = character()))

  ## Per bird, the region composition of each cluster -- so one deeply sampled bird cannot carry a
  ## label on its own.
  per_bird <- md %>%
    group_by(cluster, bird = .data[[individual_col]], region = .data[[region_col]]) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

  top <- per_bird %>%
    group_by(cluster, bird) %>%
    slice_max(frac, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(cluster, region) %>%
    summarise(n_birds = n_distinct(bird), median_frac = median(frac), .groups = "drop") %>%
    group_by(cluster) %>%
    slice_max(n_birds * median_frac, n = 1, with_ties = FALSE) %>%
    ungroup()

  batches <- if (batch_col %in% names(md)) {
    md %>%
      group_by(cluster) %>%
      summarise(n_batches = n_distinct(.data[[batch_col]]), .groups = "drop")
  } else {
    tibble(cluster = unique(md$cluster), n_batches = NA_integer_)
  }

  top %>%
    left_join(batches, by = "cluster") %>%
    mutate(
      confounded = !is.na(n_batches) & n_batches == 1,
      supported = n_birds >= min_birds & median_frac >= frac_min,
      region = if_else(supported & !confounded, region, NA_character_)
    ) %>%
    dplyr::select(cluster, region, n_birds, median_frac, n_batches, confounded, supported)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
