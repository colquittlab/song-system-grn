# Cluster filtering and naming

## Production

`hybrid_division_naming.qmd` — **the** cluster filtering and naming notebook. It is the single
place where the cluster exclusions and the published cluster labels are defined:

- input: `snrna/clustering/snrna_seurat_cellbender_preprocess/obj_clustered.qs2` (read-only)
  and `snrna/integration/composite_gg_adult/composite_calls_summary.csv`
- filtering: `excluded_clusters` (confirmed artifacts)
- naming: `celltype_hybrid` — the three excitatory divisions that hold across birds
  (`DACH2`/`SATB2`/`CACNA1H`), song-region clusters hand-labeled from anatomy/projection evidence,
  everything else a bare numbered member of its division
- output: `hybrid_division_naming/obj_hybrid_labels.qs2` and `hybrid_label_lookup.csv`

Consumers of `obj_hybrid_labels.qs2`:

- `snrna/trees/celltypes_hclust_all_hybrid.qmd`
- `snrna/trees/celltypes_hclust_glut_robustness_hybrid.qmd`
- `snrna/reduction_viz/combined_all_umap_hybrid.R`
- `snrna/deg/dotplot_manual_markers_hybrid.qmd`
- `snrna/integration/colquitt2021_label_transfer.qmd`

Change a label or an exclusion in `hybrid_division_naming.qmd`, re-render it, then re-render those.

## Not production

`exploratory/` — the naming work that led to the hybrid scheme. See `exploratory/README.md`.
Nothing outside that directory reads any of its outputs; they are cited in prose only.
