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
Results go to `<results>/config3/perturbation/` (regressors.pkl, perturb_<TF>_sv<0|10>_<iter>.parquet, PDFs); they are not tracked.

Expression scale: simulations run on the raw counts SCENIC+ itself used, so the OE value of 10 means 10 counts. Earlier runs whose `.raw`
held log-normalized values (the old double-normalization, see `ra-arco-hvc-nc/anndata_rna/make_adata_doublenorm.py`) were on a different scale; do not compare OE magnitudes across them.
