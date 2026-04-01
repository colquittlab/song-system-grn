# Song System GRN — Analysis Code

Code accompanying the manuscript:

> **Gene regulatory co-option drives birdsong neural circuit specialization**
> *(Colquitt lab — manuscript in preparation)*

---

## Data availability

Raw sequencing data are available at GEO:

| Dataset | Accession |
|---------|-----------|
| snRNA-seq | [GSE316539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316539) |
| Multiome | [GSE316538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316538) |

---

## Repository structure

### `snrna/` — snRNA-seq processing
Single-nucleus RNA-seq (CellRanger + CellBender) clustering and downstream analysis of songbird brain regions (RA, ARCO, HVC, NC, LMAN, NR).

| Directory/File | Description |
|----------------|-------------|
| `preprocessing/` | CellRanger preprocessing |
| `clustering/` | Seurat clustering with CellBender ambient RNA correction |
| `deg/` | Differential gene expression (all cell types, astrocytes, song vs. surround glutamatergic neurons) |
| `trees/` | Hierarchical clustering dendrograms (all cell types; glutamatergic) |
| `reduction_viz/` | UMAP visualizations |
| `integration/` | Integration with Zaremba et al. chicken snRNA-seq dataset |

### `multiome/` — Single-nucleus multiome (RNA + ATAC)
Joint RNA and chromatin accessibility profiling.

| Directory | Description |
|-----------|-------------|
| `seurat/` | Seurat preprocessing, clustering, DEG, and UMAP visualization |
| `archr/` | ArchR ATAC processing: peak calling, differential accessibility (DARs), ChromVAR motif enrichment, HOMER motif analysis, DREME, peak-to-gene links, positive regulator identification |
| `snapatac/` | SnapATAC2 preprocessing, peak calling, and bigwig export |

### `grn/` — Gene regulatory network inference (SCENIC+)
SCENIC+ GRN inference pipeline and downstream R analysis.

| Directory | Description |
|-----------|-------------|
| `cistarget/` | Scripts to build custom cisTarget motif-to-region databases |
| `ra-arco-hvc-nc/` | All projection neurons — SCENIC+ pipeline (config3) |
| `ra-arco-hvc-nc_glut/` | Glutamatergic projection neurons — SCENIC+ pipeline (config15) |
| `ra-arco-hvc-nc_gaba/` | GABAergic neurons — SCENIC+ pipeline (config13) |
| `ra-arco-hvc-nc_astro-oligo/` | Astrocytes + oligodendrocytes — SCENIC+ pipeline (config1) |
| `glut-ra_glut-arco-1/` | RA vs ARCO glutamatergic comparison — SCENIC+ pipeline (config18) |
| `analysis/` | R/Quarto downstream analysis of SCENIC+ eRegulon outputs |

Each SCENIC+ directory contains:
- `anndata_rna/` — AnnData preparation notebooks
- `pycisTopic/` — pycisTopic topic modeling
- `cisTarget/` — region-specific cisTarget database creation
- `scenicplus/` — Snakemake workflow + config files

### `comparative_genomics/` — Evolutionary conservation
Comparative genomics analysis across bird and mammal species.

| Directory | Description |
|-----------|-------------|
| `cactus/` | Nextflow pipelines for whole-genome alignment (Cactus), MAF extraction/filtering, phyloFit neutral model estimation, phyloP conservation scoring, phastCons; configs for avian363 and mamAvi605 alignments |
| `deeptools/` | deepTools conservation analysis at SnapATAC2 peaks; EMX2-specific analysis |

### `overexpression/` — Chicken overexpression experiments
Analysis of MAFB + EMX2 co-overexpression in chicken projection neurons.

### `scripts/` — Shared utilities
Python and shell scripts for HOMER motif annotation and MEME format conversion.

---

## Reproducing the environment

### Python

Three conda environments are used, with exported specs in `envs/`. Requires [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), [mamba](https://mamba.readthedocs.io/), or [conda](https://docs.conda.io/).

```bash
micromamba env create -f envs/scenicplus.yaml   # GRN notebooks (pycisTopic, SCENIC+)
micromamba env create -f envs/snapatac.yaml     # SnapATAC2 preprocessing
micromamba env create -f envs/scanorama.yaml    # Scanorama integration (via reticulate)
```

The Jupyter kernel display names used in the notebooks map to these environments:

| Kernel (`display_name`) | Conda environment |
|------------------------|-------------------|
| `scenicplus3` | `scenicplus` |
| `snapatac` | `snapatac` |

After creating the environments, register the kernels with Jupyter if needed:

```bash
micromamba run -n scenicplus python -m ipykernel install --user --name scenicplus --display-name "scenicplus3"
micromamba run -n snapatac   python -m ipykernel install --user --name snapatac   --display-name "snapatac"
```

### R

R package versions are captured in `renv.lock` (R 4.4.2, 563 packages). To restore:

```r
install.packages("renv")
renv::restore()
```

Some packages require additional setup before `renv::restore()`:

- **Bioconductor packages** (e.g., BiocGenerics, GenomicRanges): configure the Bioconductor repository first:
  ```r
  install.packages("BiocManager")
  BiocManager::install()  # sets up repos
  ```
- **ArchR**: installed from GitHub (`GreenleafLab/ArchR`); renv will handle this automatically but requires `remotes`.
- **Custom packages** (e.g., `TxDb.LstriataDomestica.lonStrDom2.merge3p.ucsc.assembled`): contact the Colquitt lab for access.

### Reference data

Edit `config/paths.R` (R/Quarto scripts) and `config/paths.sh` (shell scripts) to point to your local copies of the reference files. See `data/README.md` for a full list with download sources.

---

## Software

| Tool | Version | Reference |
|------|---------|-----------|
| CellRanger | 7.x | 10x Genomics |
| CellBender | 0.3.x | Fleming et al. 2023 |
| Seurat | 5.x | Hao et al. 2024 |
| ArchR | 1.0.x | Granja et al. 2021 |
| SnapATAC2 | 2.x | Zhang et al. 2024 |
| SCENIC+ | 1.0.x | Bravo González-Blas et al. 2023 |
| Cactus | 2.6.x | Armstrong et al. 2020 |
| deepTools | 3.x | Ramírez et al. 2016 |
| HOMER | 4.x | Heinz et al. 2010 |
| R | 4.3+ | R Core Team |
| Python | 3.10+ | |
