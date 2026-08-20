# `xenium/` — Xenium in situ spatial transcriptomics

Spatial localization of the snRNA-seq-defined cell types (`snrna/clustering/`)
in the adult Bengalese finch brain.

## Dataset

Run `Brainard_KL_052905` (2025-05-29), Xenium Ranger 3.3.0.1, fresh-frozen,
custom 425-gene panel (design ID `Z2XBPV`, "otherBF Brain_425g").

**Ten sections — two birds × five coronal levels**, ~1.49M cells:

| Bird / cassette | Slide | Regions | Cells per region |
|---|---|---|---|
| OR52YW26 | 0028723 | 1_4, 1_7, 2_2, 2_4, 2_7 | 181k, 179k, 165k, 145k, 134k |
| OR69PU4 | 0028750 | 1_4, 1_7, 2_2, 2_4, 2_7 | 165k, 185k, 139k, 133k, 167k |

Median 412 transcripts / 135 genes per cell; 83% of transcripts assigned to
cells; adjusted negative-control probe rate 0.36%.

Paths are configured in `config/paths.R` (`XENIUM_DIR`, `XENIUM_RAW_DIR`,
`XENIUM_PANEL_CSV`, `XENIUM_HE_DIR`, `XENIUM_SELECTIONS_DIR`).

### Dataset caveats

- **Segmentation is nucleus expansion only.** Every section reports
  `segmented_cell_nuc_expansion_frac = 1.0` — no boundary stain was used. Expect
  transcript spillover between neighbouring cells; prefer mixture-aware
  assignment (RCTD) over hard nearest-reference assignment, and retain
  `transcripts.parquet` in case resegmentation is needed.
- **Region labels in the collaborator-supplied objects are unreliable.** The
  section comments in `Xenium Codes/Xenium.R` flag mismatches (e.g.
  `##### HVC 1_7 (Actually 2_2) ######`). Take `region_name` from each region's
  `experiment.xenium`, and rebuild objects from the raw outputs rather than
  inheriting `Xenium Objects/*.RData`.
- **EMX2 and MAFB are not on the panel.** They can only be read indirectly,
  through their SCENIC+ target modules.

## Contents

| Directory | Description |
|-----------|-------------|
| `label_transfer/` | Panel resolution benchmark; transfer of snRNA-seq cell type labels onto the Xenium cells |
| `preprocessing/` | Per-region loading, QC, normalization, and merging of the ten sections |
| `spatial_analysis/` | Song nucleus vs. surround composition, boundary sharpness, niche analysis, eRegulon spatial scoring |

### `preprocessing/xenium_preprocess.qmd`

Builds the analysis object: ten sections loaded from the Xenium Ranger outputs
(centroids only — molecule coordinates and segmentation polygons stay on disk),
per-cell covariates joined from `cells.csv.gz`, QC filtered at ≥20 transcripts
and ≥10 genes, merged into one object with one field of view per section, and
log-normalized to the median transcript count — the same normalization the panel
benchmark applies to the reference, so the two match at label transfer.

**Result:** 1,586,249 of 1,594,009 cells retained (99.2–99.8% per section). Mean
control fraction 0.043–0.059%, so background is negligible. Harmony integration
over `section_id` gives 24 QC clusters with sections and birds fully intermixed —
no cluster is batch-driven — and those clusters form coherent anatomical domains
when plotted back into tissue, confirming the merge and per-section coordinate
frames. The `SLC17A6` / `GAD1` complementarity across the pallium–striatum
boundary reproduces on every section.

The embedding runs on a stratified subsample (10k cells per section). Seurat's
`FindNeighbors` and `FindClusters` are single-threaded with no thread argument
and take >2 h on the full 1.6M cells, and the clusters they produce are a sanity
check that transferred labels replace; `ScaleData`/PCA/Harmony use `future` and
OpenBLAS threading, and `RunUMAP` passes `n_threads` through to uwot. The
marker-gene spatial panels use **all** cells, so the anatomical check does not
depend on the subsample.

Outputs — objects in `XENIUM_PROC_DIR/xenium_preprocess/`:

| File | Description |
|------|-------------|
| `obj_xenium.qs2` | Full object: counts, normalized data, coordinates, QC clusters on the subsampled cells (839 MB) |
| `obj_xenium_qc-subsample.qs2` | 100k-cell subsample with PCA / Harmony / UMAP reductions (118 MB) |

Tables and figures in `xenium/preprocessing/xenium_preprocess/`:

| File | Description |
|------|-------------|
| `sections.csv` | Run metadata per section, from each `experiment.xenium` |
| `qc_summary_prefilter.csv`, `qc_retention.csv` | Per-section QC and filter retention |
| `qc_distributions.pdf` | Counts, genes, cell area, control fraction per section |
| `umap_qc.png`, `qc_cluster_composition.*` | QC embedding and per-section cluster composition |
| `spatial_cluster_qc_<section>.png` | QC clusters in tissue, per section |
| `spatial_markers_<section>.png` | Marker genes on all cells, per section |
| `cell_metadata.csv.gz` | Per-cell metadata for the full object |

### `label_transfer/panel_resolution_benchmark.qmd`

Establishes the label granularity the panel can support, and therefore what
every downstream analysis is allowed to assume. Run this first.

The snRNA-seq reference is subset to the 425 panel genes and binomially thinned
to the transcript depth observed in the Xenium run, then held-out cells are
classified under 5-fold cross-validation at three levels of the taxonomy
(`cluster` → `subclass` → `class`) with two independent classifiers:
multinomial naive Bayes (measures the panel's information content) and Seurat
transfer anchors (measures what the applied workflow recovers). A recommended
label set is then derived by iteratively merging the worst-performing label with
whichever label carries the most confusion mass with it, until every label
clears an F1 threshold — a merge driven by the confusion structure rather than
by the naming hierarchy.

**Result:** all 425 panel genes are present in the reference feature space, and
the taxonomy survives the panel nearly intact. Macro-F1 at the full 50-cluster
resolution is 0.85 (naive Bayes) / 0.80 (transfer anchors); 0.95 at class level
under both. Only one merge is required (`Glut-NC-1` + `Glut-Nido-1` →
`Glut-NC/Glut-Nido`), giving a recommended set of **49 labels** at macro-F1 0.86
/ 0.81. Song-nucleus types are among the strongest: Glut-RA 0.98, Glut-HVC-2
0.96, Glut-LMAN-2 0.96, Glut-HVC-1 0.95, Glut-LMAN-1 0.89 (naive Bayes F1).
Weakest survivors are Glut-Arco-3 (0.61) and Glut-Arco-1 (0.69), which confuse
with each other, and GABA-Pre (0.69). GABA-Pre and Glut-Arco-3 fall below
threshold under transfer anchors (0.29, 0.48) while passing under naive Bayes —
i.e. the panel carries the information but anchor-based transfer loses it, an
argument for a mixture-aware method (RCTD) on the real data.

Outputs (`xenium/label_transfer/panel_resolution_benchmark/`):

| File | Description |
|------|-------------|
| `panel_genes_in_reference.csv` | Panel gene presence in the reference feature space |
| `cv_predictions.qs2` | Per-cell cross-validated predictions, both classifiers, all levels |
| `per_label_metrics.csv` | Precision / recall / F1 per label, per level, per classifier |
| `summary_metrics.csv` | Macro-F1 and balanced accuracy per level and classifier |
| `recommended_label_set.csv` | Per-cluster verdict and the recommended transfer target |
| `merge_log.csv` | Which labels were merged, at what F1, and with what |
| `per_label_metrics_recommended.csv`, `summary_metrics_recommended.csv` | Metrics for the recommended label set |
| `confusion_*.pdf`, `f1_by_method.pdf` | Confusion structure and method comparison |
