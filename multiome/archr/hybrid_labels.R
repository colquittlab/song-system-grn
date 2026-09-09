# multiome/archr/hybrid_labels.R
#
# Attach the hybrid-scheme cluster labels to an ArchR project, in memory.
#
#   source(here::here("config/paths.R"))
#   source(here::here("multiome/archr/hybrid_labels.R"))
#
#   proj = loadArchRProject(ARCHR_PROJ_DIR)
#   proj = add_cluster_hybrid(proj)      # adds proj$cluster_hybrid, drops unlabeled cells
#
# The labels come from multiome/seurat/label_transfer_hybrid.qmd, which transferred
# the hybrid snRNA labels (snrna/naming/hybrid_division_naming.qmd, `celltype_hybrid`)
# onto the multiome object and curated them to one label per multiome cluster by
# majority vote. That curation is a per-cluster map, so it is applied here from the
# tracked CSV rather than by loading the 2 GB `obj_clustered_hybrid.qs2` --
# `verify_cluster_hybrid_map()` below checks the two agree cell-by-cell when the
# object is available.
#
# Nothing here writes to the ArchR project on disk. `add_cluster_hybrid()` modifies
# only the in-memory `cellColData`; call `saveArchRProject()` yourself if you mean to
# persist it (no script in this directory does).

HYBRID_LABEL_DIR <- here::here("multiome/seurat", "label_transfer_hybrid")
HYBRID_VOTE_CSV  <- file.path(HYBRID_LABEL_DIR, "cluster_hybrid_majority_vote.csv")
HYBRID_OBJ_QS2   <- file.path(HYBRID_LABEL_DIR, "obj_clustered_hybrid.qs2")

## Curated multiome cluster -> hybrid label. Mostly a 1:1 rename; the two real
## regroupings are Glut-Arco-2/6/7 -> Glut-CACNA1H-2 and Astro-1/2 -> Astro.
## `Glut-Nido-3` is absent by design -- confirmed artefactual in the label-transfer
## notebook (scattered across the NC-1/NC-2/NC-4 territory rather than its own) and
## excluded there, so it has no hybrid label and is dropped below.
hybrid_cluster_map <- function() {
  vote <- readr::read_csv(HYBRID_VOTE_CSV, show_col_types = FALSE)
  setNames(vote$cluster_hybrid, vote$cluster)
}

## Order used for every hybrid-scheme figure: the multiome's own subset of the
## division/family ordering in snrna/trees/celltypes_hclust_all_hybrid.qmd.
hybrid_ct_order <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx",
                     "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
                     "Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2",
                     "Glut-Im", "Glut-NB", "Glut-NSC",
                     "GABA-LGE-1", "GABA-LGE-2",
                     "GABA-MGE-SST-1", "GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2",
                     "GABA-MGE-LAMP5", "GABA-MGE-LHX8",
                     "GABA-CGE", "GABA-Im",
                     "Astro", "Epen", "Oligo", "OPC", "Micro", "Endo")

## The focal song-vs-surround contrasts, carried over from the pre-hybrid scripts
## (Glut-RA/Glut-Arco-1, Glut-HVC-1/Glut-NC-1, Glut-HVC-2/Glut-NC-4) under the new
## names. Same cells on both sides -- only the labels changed.
hybrid_pairs_glut <- list(ra    = list("Glut-CACNA1H-RA", "Glut-CACNA1H-1"),
                          hvc_1 = list("Glut-DACH2-HVCra", "Glut-DACH2-1"),
                          hvc_2 = list("Glut-DACH2-HVCx", "Glut-DACH2-4"))

## GABA-4-1 vs GABA-1-1 and GABA-2-1 vs GABA-1-1, renamed.
hybrid_pairs_gaba <- list(pvalb2_lge1 = list("GABA-MGE-PVALB-2", "GABA-LGE-1"),
                          sst1_lge1   = list("GABA-MGE-SST-1", "GABA-LGE-1"))

## Semantic groupings the pre-hybrid scripts built with grepl() over the old names
## ("Arco|RA", "HVC|NC|Nido", "Pre"). Those patterns don't survive the rename, so the
## groups are spelled out instead of re-derived from a pattern.
hybrid_groups <- list(
  glut       = grep("^Glut-", hybrid_ct_order, value = TRUE),
  gaba       = grep("^GABA-", hybrid_ct_order, value = TRUE),
  ## Neuronal precursors, both divisions (was Glut-Pre-1/2/3 + GABA-Pre).
  precursor  = c("Glut-Im", "Glut-NB", "Glut-NSC", "GABA-Im"),
  ## The excitatory side of the Glut-vs-GABA contrast in peaks_gaba-glut.qmd. That
  ## script excluded Glut-Pre-1 and Glut-Pre-2 but kept Glut-Pre-3; under the hybrid
  ## names that is everything Glut except Glut-NSC and Glut-NB, with Glut-Im kept.
  glut_broadclass = c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx",
                      "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
                      "Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2",
                      "Glut-Im"),
  ## RA/arcopallial family (was Glut-RA + Glut-Arco-*).
  arco_glut  = c("Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2"),
  ## HVC/nidopallial family (was Glut-HVC-* + Glut-NC-* + Glut-Nido-3).
  nido_glut  = c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx",
                 "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4"),
  neuroglia  = c("Astro", "Oligo", "OPC")
)

## Dissection each focal cluster belongs to, keyed for `position_colors2`. The
## pre-hybrid scripts read this straight out of the label ("Glut-HVC-1" -> "hvc" by
## splitting on "-"), which the hybrid names no longer support -- "Glut-DACH2-HVCra"
## splits to "DACH2". Spelled out for the six song/surround clusters that use it.
hybrid_position <- c("Glut-DACH2-HVCra" = "hvc",
                     "Glut-DACH2-HVCx"  = "hvc",
                     "Glut-DACH2-1"     = "nc",
                     "Glut-DACH2-4"     = "nc",
                     "Glut-CACNA1H-RA"  = "ra",
                     "Glut-CACNA1H-1"   = "arco")

## Song nuclei among the focal clusters (the rest are their surrounds).
hybrid_song_clusters <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-CACNA1H-RA")

#' Add `cluster_hybrid` to an ArchR project's cellColData, in memory.
#'
#' @param proj ArchRProject carrying the curated `cluster` column.
#' @param drop_unlabeled drop cells whose `cluster` has no hybrid label
#'   (i.e. Glut-Nido-3). TRUE matches the label-transfer notebook, which excludes
#'   them from its saved object.
#' @param as_factor store `cluster_hybrid` as a factor levelled by `hybrid_ct_order`.
#'   getMarkerFeatures() wants a character column, so the default is FALSE.
add_cluster_hybrid <- function(proj, drop_unlabeled = TRUE, as_factor = FALSE) {
  map <- hybrid_cluster_map()
  clusters <- as.character(proj$cluster)

  unlabeled <- setdiff(unique(clusters), names(map))
  if (length(unlabeled) > 0 && !identical(unlabeled, "Glut-Nido-3")) {
    ## Anything beyond the known-artefactual cluster means the ArchR project and the
    ## label-transfer notebook have drifted apart -- don't silently drop those cells.
    stop("Cluster(s) with no hybrid label, and not the expected Glut-Nido-3: ",
         paste(unlabeled, collapse = ", "))
  }

  hybrid <- unname(map[clusters])
  proj$cluster_hybrid <- hybrid

  if (drop_unlabeled && anyNA(hybrid)) {
    keep <- getCellNames(proj)[!is.na(hybrid)]
    message("Dropping ", sum(is.na(hybrid)), " cell(s) with no hybrid label (",
            paste(unlabeled, collapse = ", "), "); ", length(keep), " remain.")
    proj <- subsetCells(proj, keep)

    ## subsetCells() drops the cells from cellColData but leaves the stored
    ## embeddings and reduced dimensions at their original size, so plotEmbedding()
    ## then refuses with "Not all cells in embedding are present in ArchRProject!".
    ## Carry the subset through to them as well.
    for (nm in names(proj@embeddings)) {
      d <- proj@embeddings[[nm]]$df
      proj@embeddings[[nm]]$df <- d[rownames(d) %in% keep, , drop = FALSE]
    }
    for (nm in names(proj@reducedDims)) {
      m <- proj@reducedDims[[nm]]$matDR
      if (!is.null(m)) {
        proj@reducedDims[[nm]]$matDR <- m[rownames(m) %in% keep, , drop = FALSE]
      }
    }

    ## The impute weights are cell-by-cell matrices held in files inside the shared
    ## project, so they cannot be subset without rewriting them there. Drop them
    ## instead: plotEmbedding() then plots unimputed values rather than failing on
    ## the size mismatch. Call addImputeWeights(proj) on the subset project if
    ## smoothed values are wanted.
    if (length(proj@imputeWeights) > 0) {
      message("Dropping the project's impute weights -- they cover the full cell set. ",
              "Call addImputeWeights(proj) if you want imputed values.")
      proj@imputeWeights <- S4Vectors::SimpleList()
    }
  }

  present <- setdiff(unique(as.character(proj$cluster_hybrid)), hybrid_ct_order)
  if (length(present) > 0) {
    stop("Hybrid label(s) missing from hybrid_ct_order: ", paste(present, collapse = ", "))
  }
  if (as_factor) {
    proj$cluster_hybrid <- factor(as.character(proj$cluster_hybrid),
                                  levels = intersect(hybrid_ct_order,
                                                     unique(as.character(proj$cluster_hybrid))))
  }
  proj
}

#' Check the per-cluster map against the label-transfer object's per-cell labels.
#'
#' The map is only a shortcut if every cell in `obj_clustered_hybrid.qs2` carries the
#' label the map would give it. Run once after the label transfer is re-run; needs
#' qs2 and ~2 GB of the object in memory, which is why it is not called on import.
verify_cluster_hybrid_map <- function() {
  stopifnot(file.exists(HYBRID_OBJ_QS2))
  obj <- qs2::qs_read(HYBRID_OBJ_QS2, nthreads = 8)
  md <- obj@meta.data
  map <- hybrid_cluster_map()
  expected <- unname(map[as.character(md$cluster)])
  ok <- identical(expected, as.character(md$cluster_hybrid))
  if (!ok) {
    mism <- which(expected != as.character(md$cluster_hybrid))
    stop(length(mism), " cell(s) disagree with the per-cluster map, e.g. ",
         paste(utils::head(rownames(md)[mism], 3), collapse = ", "))
  }
  message("cluster_hybrid map reproduces all ", nrow(md), " per-cell labels in ",
          basename(HYBRID_OBJ_QS2), ".")
  invisible(TRUE)
}
