# ra-arco-hvc-nc_gaba_hybrid

SCENIC+ on **GABAergic neurons**, hybrid labels (`cluster_hybrid`). Derived from `ra-arco-hvc-nc_glut_hybrid` (see that README) by
`../make_hybrid_subset_dirs.py`, changing only what the earlier run (ra-arco-hvc-nc_gaba (config13)) specified.

| | earlier run | this directory |
|---|---|---|
| cells | older labels | all `GABA-*` clusters minus `GABA-Im` (the renamed `GABA-Pre`): 8 clusters, 3,275 RNA cells |
| `.raw` | older export (not re-checked here) | raw counts |
| LDA topics | 20 (hand-picked from 2-40) | 20 (only this model trained) |
| DARs (every cluster) | adj. p 0.05, abs(log2FC) >= log2(1.5) | same; no song-pair DARs (the earlier run had none) |
| motif stage | motif_similarity_fdr 0.01, dem_adj_pval_thr 0.05, dem_log2fc_thr 1.0, ctx_auc_threshold 0.01, ctx_nes_threshold 2.0, ctx_rank_threshold 0.1 | **config1 = this**, config2 = the all-cell strict thresholds; both with the hybrid cisTarget database |

## Steps (same as the glut run)

1. `anndata_rna/make_gaba_adata.py` -> `anndata_rna/adata.h5ad` (`.raw` = raw counts).
2. `pycisTopic/pycistopic_gaba.py` (script) or `pycisTopic/pycisTopic.ipynb` (the same code, interactive): subset the all-cell hybrid cisTopic object to
   the cells with ATAC, train the 20-topic model, write topic sets and DARs. Output `cistopic_obj_gaba.pkl` and `region_sets/`.
3. `scenicplus/make_configs.py` -> config1 and config2, **both on the hybrid cisTarget database** (`ra-arco-hvc-nc_hybrid`, built on the cisTopic object's own
   consensus regions; there is no config on the older shared database here). Config2 reuses config1's accessibility+expression object, search space and
   region-to-gene fit: run `scenicplus/submit_all.sh 1` on prism, then `scenicplus/prepare_shared_upstream.sh`, then `scenicplus/submit_all.sh 2`.
   **Configs 3 and 4** test the region-to-gene search space: config2 with up to 250 kb (config3) or 500 kb (config4) upstream and downstream of each
   gene, instead of the default 150 kb. They reuse only config1's accessibility+expression object (`shared/outs/ACC_GEX.h5mu`) and compute their own
   search space, region-to-gene fit and everything after it (12 h job limits, the medium partition's maximum). With `shared/outs` in place on prism:
   `scenicplus/submit_all.sh 3 4`. Same design as `glut-ra_glut-c1h-1_allruns_hybrid` configs 3-4; compare 2 -> 3 -> 4.
   Helpers (`submit_all.sh`, `link_results_to_hdd.sh`, `prepare_shared_upstream.sh`) are the glut run's with this root.

Data live under `/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/.../ra-arco-hvc-nc_gaba_hybrid/`, outside git. Analysis templates are not made yet (they need results).
