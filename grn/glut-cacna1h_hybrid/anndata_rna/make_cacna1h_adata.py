#!/usr/bin/env python3
"""all Glut-CACNA1H subtypes (Glut-CACNA1H-RA, -1 and -2): SCENIC+ expression input, hybrid labels (a new subset, no earlier run; modelled on glut-ra_glut-arco-1).

Same steps as ../../ra-arco-hvc-nc/anndata_rna/anndata.ipynb (which reads data/adata.h5ad from export_h5ad_hybrid.R), with the cell filter of the
earlier run: all clusters containing "CACNA1H" (Glut-CACNA1H-RA, Glut-CACNA1H-1, Glut-CACNA1H-2; older labels Glut-RA, Glut-Arco-1, Glut-Arco-2 and relatives) from ALL
runs. There is no earlier run of this subset to copy a run3 exclusion from (the RA vs Arco-1 run excluded run3; the two RA vs C1H-1 directories here
show both options), so every run is kept.

    adata.h5ad   .raw = raw integer counts (what SCENIC+ reads); X = log-normalized, HVG, scaled (display only)

    python make_cacna1h_adata.py
"""
import os

import anndata as ad
import numpy as np
import scanpy as sc

B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/")
SRC = B + "ra-arco-hvc-nc_hybrid/data/adata.h5ad"          # X = counts, cluster = cluster_hybrid
OUT_DIR = B + "glut-cacna1h_hybrid/anndata_rna/"
CLUSTERS = ["Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2"]

adata = sc.read_h5ad(SRC)
assert np.all(adata.X.data == np.round(adata.X.data)), "source X is not integer counts"

cl = adata.obs["cluster"].astype(str)
keep = cl.str.contains("CACNA1H")   # all runs
adata = adata[keep.values].copy()
assert set(adata.obs["cluster"].astype(str)) == set(CLUSTERS), "unexpected CACNA1H clusters: " + str(sorted(set(adata.obs["cluster"].astype(str))))
adata.obs["cluster"] = adata.obs["cluster"].astype(str)
print(adata.obs["cluster"].value_counts().to_string())

adata.obsm["archr_umap"] = adata.obsm["X_dims30nn30mindist0.3"]

# .raw must stay raw counts; `adata.raw = adata` shares the sparse X, so copy before the in-place normalize/log1p.
adata.raw = adata.copy()
assert np.all(adata.raw.X.data == np.round(adata.raw.X.data)), ".raw is not integer counts"
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)

# barcodes -> the cisTopic object's "<barcode>-1___<sample_short>" names
obs = adata.obs.copy()
obs.index = obs.index.str.replace("-[0-9]", "-1", regex=True)
obs["barcode"] = obs.index
obs.index = obs.apply(lambda x: "___".join([x.barcode, f"{x.sample_short}"]), axis=1)
assert obs.index.is_unique
adata.obs = obs

os.makedirs(OUT_DIR, exist_ok=True)
adata.write_h5ad(OUT_DIR + "adata.h5ad")
chk = ad.read_h5ad(OUT_DIR + "adata.h5ad")
print(chk, "| raw:", chk.raw.shape, "integer:", bool(np.all(chk.raw.X.data == np.round(chk.raw.X.data))))
