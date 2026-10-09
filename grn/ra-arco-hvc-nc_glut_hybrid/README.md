# ra-arco-hvc-nc_glut_hybrid

SCENIC+ on the **glutamatergic projection neurons only**, with the **hybrid cluster labels** (`cluster_hybrid`).
It repeats the `ra-arco-hvc-nc` hybrid pipeline on the subset that `ra-arco-hvc-nc_glut` defined for the older labels.

| | `ra-arco-hvc-nc_glut` (older labels) | `ra-arco-hvc-nc` (hybrid, all cells) | this directory |
|---|---|---|---|
| labels | seurat-snrna-clustering | hybrid | hybrid |
| cells | `Glut` and not `Pre` | all 27 clusters | `Glut` minus `Glut-NSC/NB/Im` (the renamed `Glut-Pre-1/2/3`): 9 clusters, 9,159 RNA cells (9,063 with ATAC) |
| `.raw` | older export (not re-checked here) | raw counts | raw counts |
| LDA topics | 20 | 20 (also 15, 30, 40) | 20 (only this model trained) |
| motif stage | config15: fdr 0.001, dem 0.1/0.5, auc 0.0025, nes 2.0, rank 0.2 | config1: 0.01, 0.05/1.0, 0.005, 3.0, 0.2 | **config1 = the glut spec**, config2 = the all-cell strict thresholds |
| cisTarget DB | shared `ra-arco-hvc-nc_seurat-clustering` | same | same (not rebuilt; this run's regions are a subset of its columns) |

The precursor stages are listed explicitly rather than matched on `Pre`: the hybrid relabel (`snrna/naming/hybrid_division_naming.qmd`)
renamed them, so the old substring filter would silently keep them.

## Steps

1. `anndata_rna/make_glut_adata.py` -- glut subset of `ra-arco-hvc-nc_hybrid/data/adata.h5ad` (from `export_h5ad_hybrid.R`);
   `.raw` = raw counts. Writes `anndata_rna/adata.h5ad`.
2. `pycisTopic/pycistopic_glut.py` -- subset the all-cell hybrid cisTopic object to the glut cells that have ATAC, train the 20-topic LDA model,
   write topic sets and DARs for every cluster (pycisTopic, as in the all-cell run). **Song-pair DARs are the ArchR ones from the all-cell
   run, copied**, not recomputed: `DARs_song-pairs` holds exactly Glut-CACNA1H-RA vs -1, Glut-DACH2-HVCra vs -1 and Glut-DACH2-HVCx vs -4
   (`ra-arco-hvc-nc_hybrid/pycisTopic/region_sets_k40_archr/DARs_song-pairs`, from `make_song_pair_dars_archr.R`: ArchR's bias-matched
   Wilcoxon on the peak matrix). pycisTopic's `find_diff_features` ranks imputed accessibility, which calls DARs from topics the foreground
   barely uses; the ArchR test does not, and the pairs do not depend on the cell set or topic count.
3. `scenicplus/make_configs.py` -- config1 (glut spec) and config2 (strict control); renders through `ra-arco-hvc-nc/scenicplus/make_configs.py`
   so the YAML and the per-config `run_snakemake.sbatch` are the all-cell run's (only job name and paths differ). The helpers
   `submit_all.sh`, `link_results_to_hdd.sh` (store: `.../ra-arco-hvc-nc_glut_hybrid/results`) and `.gitignore` are copies of the
   all-cell ones. On prism, after copying `anndata_rna/adata.h5ad`, `pycisTopic/cistopic_obj_glut.pkl` and `pycisTopic/region_sets/`
   there: `scenicplus/submit_all.sh` (all configs) or `scenicplus/submit_all.sh 1` (config1 only).
4. `../analysis/scenicplus_hybrid_glut_template.qmd` -- the all-cell hybrid analysis restricted to the glut object, with the song-pair DEG
   sets `ra`, `hvcra`, `hvcx` in place of `gaba4`/`astro`.

Data live under `/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/.../ra-arco-hvc-nc_glut_hybrid/`, outside git.
