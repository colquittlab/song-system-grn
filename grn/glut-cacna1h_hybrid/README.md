# glut-cacna1h_hybrid

SCENIC+ on **all Glut-CACNA1H subtypes (Glut-CACNA1H-RA, -1 and -2)**, hybrid labels (`cluster_hybrid`). Derived from `ra-arco-hvc-nc_glut_hybrid` (see that README) by
`../make_hybrid_subset_dirs.py`, changing only what the earlier run (glut-ra_glut-arco-1 (config18; no earlier run of this subset: its thresholds and topic count are reused)) specified.

| | earlier run | this directory |
|---|---|---|
| cells | older labels | all three `Glut-CACNA1H` subtypes (RA 560, CACNA1H-1 1,406, CACNA1H-2 1,433) from all runs: 3,399 RNA cells. New subset, no earlier run: thresholds and topic count (10) are those of the RA vs C1H-1 run, not re-selected. |
| `.raw` | older export (not re-checked here) | raw counts |
| LDA topics | 10 (hand-picked from 2-40) | 10 (only this model trained) |
| DARs (every cluster) | adj. p 0.1, abs(log2FC) >= log2(1.5) | same; the ArchR song pair copied from the all-cell run (`Glut-CACNA1H-RA` vs `Glut-CACNA1H-1`) is kept in `pycisTopic/archr_song_pair/`, NOT in `region_sets/`: DEM tests every set against the others in its folder, so a one-set folder has an empty background and fails (`No objects to concatenate`); `DARs_all` already has the RA-up and CACNA1H-1-up sets |
| motif stage | motif_similarity_fdr 0.001, dem_adj_pval_thr 0.05, dem_log2fc_thr 1.0, ctx_auc_threshold 0.005, ctx_nes_threshold 3.0, ctx_rank_threshold 0.05; extended_annotation 'Orthology_annot Motif_similarity_annot' | **config1 = this**, config2 = the all-cell strict thresholds; both with the hybrid cisTarget database |

## Steps (same as the glut run)

1. `anndata_rna/make_cacna1h_adata.py` -> `anndata_rna/adata.h5ad` (`.raw` = raw counts).
2. `pycisTopic/pycistopic_cacna1h.py` (script) or `pycisTopic/pycisTopic.ipynb` (the same code, interactive): subset the all-cell hybrid cisTopic object to
   the cells with ATAC, train the 10-topic model, write topic sets and DARs. Output `cistopic_obj_cacna1h.pkl` and `region_sets/`.
3. `scenicplus/make_configs.py` -> config1 and config2, **both on the hybrid cisTarget database** (`ra-arco-hvc-nc_hybrid`, built on the cisTopic object's own
   consensus regions; there is no config on the older shared database here). Config2 reuses config1's accessibility+expression object, search space and
   region-to-gene fit: run `scenicplus/submit_all.sh 1` on prism, then `scenicplus/prepare_shared_upstream.sh`, then `scenicplus/submit_all.sh 2`.
   Helpers (`submit_all.sh`, `link_results_to_hdd.sh`, `prepare_shared_upstream.sh`) are the glut run's with this root.

Data live under `/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/.../glut-cacna1h_hybrid/`, outside git. Analysis templates are not made yet (they need results).
