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

**Flow-field figures.** `PERTURB_ARROW_STYLE` = `stream` (default, SCENIC+'s own streamplot, written to `pca_lines_*.pdf`), `quiver` (individual arrows at the grid
points, `pca_arrows_*.pdf`; `PERTURB_ARROW_MIN_FRAC`, default 0.02, keeps arrows at least that fraction of the longest) or `both`. SCENIC+'s streamplot drops arrows
shorter than 15% of the longest, so a sparse field (a TF expressed in few cells, e.g. MAFB at 7% of cells) can come out empty although the field exists; the quiver shows it.
Each cell's arrow is built from a kNN of 20% of all cells in the 2-D Harmony-PCA embedding, randomly subsampled with unseeded `np.random.choice`: set `PERTURB_SEED` to make it
reproducible. `PERTURB_FIG_SIZE` (default 1.8) is the figure side in inches and `PERTURB_POINT_SIZE` (default 1.2) and `PERTURB_POINT_ALPHA` (default 0.7) the cell scatter size and opacity (semi-transparent, so arrows show through; quiver arrows above 10% of the longest get a thin white outline); cluster labels are 5 pt and titles 6 pt Arial (embedded TrueType); cluster colours are tab10 without the gray (tab20 beyond 9 clusters). With `PERTURB_SEED` set the quiver's arrow field is cached next to the figures (`arrow_field_<kd|oe>_<TF>_seed<N>.pkl`, ~10 min to compute), so restyling does not recompute it.
`PERTURB_REDO_PERTURBATIONS=0` redraws plots from cached simulations; `PERTURB_OE_TFS=none|auto|<TFs>` and `PERTURB_OE_SET_VALUE` control overexpression
(its value is on the scale of the expression matrix); `PERTURB_EXPRESSION_SCALE=lognorm` is for a matrix built with `assemble_perturbation_h5mu.py --lognorm`.

Both calls to the child script use the notebook's own interpreter (`sys.executable`) and a non-zero exit status stops the run; the copy this was derived from used a bare `python` and ignored failures, so a wrong environment produced "Done" with no results.

Expression scale: simulations run on the raw counts SCENIC+ itself used, so the OE value of 10 means 10 counts. Earlier runs whose `.raw`
held log-normalized values (the old double-normalization, see `ra-arco-hvc-nc/anndata_rna/make_adata_doublenorm.py`) were on a different scale; do not compare OE magnitudes across them.
