# SCENIC+ in-silico perturbation (glut, hybrid labels)

Knockout and overexpression simulations from a config's eRegulons, with the glut run's `perturbation_parallel_template.ipynb`
(analysis cells unchanged; see the notebook's first cell for what was parameterized).

| File | What |
|---|---|
| `perturbation_parallel_template.ipynb` | Batch notebook: per config, trains GBM regressors for the top-5,000 eRegulon target genes, simulates KO (set value 0) and OE (10) of `KO_TFS` over 5 iterations, plots predicted log2FC and the Harmony PCA flow field. Outputs stripped. |
| `perturbation_parallel.py` | Per-TF parallel simulation, called by the notebook (copy of the script that lived outside the repo). |
| `assemble_perturbation_h5mu.py` | Builds `scplusmdata_perturb.h5mu` for a config run with `build_scplus_mudata: false` (configs 3/4 here): the real `scRNA_counts` (raw integer counts of `anndata_rna/adata.h5ad` `.raw`, as in a full `scplusmdata.h5mu`) + that config's AUC matrices and eRegulon tables; placeholder ATAC. |
| `run_notebook.py` | Executes the notebook's code cells without nbconvert/papermill. |

```
python assemble_perturbation_h5mu.py <results>/config3/outs
PERTURB_INPUT_DIR=<results> PERTURB_CONFIGS=config3 PERTURB_H5MU=scplusmdata_perturb.h5mu PERTURB_N_JOBS=40 python run_notebook.py
```
Results go to `<results>/config3/perturbation/` (`PERTURB_OUT_NAME` changes the folder name) (regressors.pkl, perturb_<TF>_sv<0|10>_<iter>.parquet, PDFs); they are not tracked.

`PERTURB_KO_TFS=MAFB` restricts the TFs; `PERTURB_FORCE_TARGETS_OF=MAFB` gives every eRegulon target of that TF a regressor, on top of the top-5,000 genes (a TF with a small eRegulon otherwise reaches only the fraction of its targets that rank in the top 5,000). An existing `regressors.pkl` is topped up with just the missing genes, so a variant can be seeded by copying it.

Both calls to the child script use the notebook's own interpreter (`sys.executable`) and a non-zero exit status stops the run; the copy this was derived from used a bare `python` and ignored failures, so a wrong environment produced "Done" with no results.

Expression scale: simulations run on the raw counts SCENIC+ itself used, so the OE value of 10 means 10 counts. Earlier runs whose `.raw`
held log-normalized values (the old double-normalization, see `ra-arco-hvc-nc/anndata_rna/make_adata_doublenorm.py`) were on a different scale; do not compare OE magnitudes across them.
