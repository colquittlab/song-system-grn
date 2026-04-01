# envs/

Conda environment specifications for reproducing the Python analyses.

| File | Used by | Kernel / environment |
|------|---------|----------------------|
| `scenicplus.yaml` | All GRN notebooks (`grn/*/anndata_rna/`, `grn/*/pycisTopic/`) and SCENIC+ Snakemake workflows | `scenicplus3` kernel → `scenicplus` conda env |
| `snapatac.yaml` | SnapATAC2 notebooks (`multiome/snapatac/`) | `snapatac` kernel → `snapatac` conda env |
| `scanorama.yaml` | Seurat Scanorama integration (called via `reticulate` from `multiome/seurat/seurat_preprocess.qmd` and `snrna/clustering/`) | `scanorama` conda env |

## Restoring an environment

```bash
micromamba env create -f envs/scenicplus.yaml
micromamba env create -f envs/snapatac.yaml
micromamba env create -f envs/scanorama.yaml
```

Or with conda/mamba:

```bash
conda env create -f envs/scenicplus.yaml
```

## R environment

R package versions are captured in `renv.lock` (R 4.4.2, 563 packages).
To restore:

```r
install.packages("renv")
renv::restore()
```

Note: some packages are installed from GitHub (e.g., ArchR) and Bioconductor.
The `renv.lock` records the source repository and version for each package.
