#!/usr/bin/env python3
"""Glutamatergic-only SCENIC+ expression input, hybrid labels (the hybrid counterpart of ra-arco-hvc-nc_glut/anndata_rna/anndata.ipynb).

Same steps as ../../ra-arco-hvc-nc/anndata_rna/anndata.ipynb (which reads data/adata.h5ad from export_h5ad_hybrid.R), with the
glut cell filter of the earlier glut run: clusters containing "Glut", minus the precursor stages. That run excluded them with
`~cluster.str.contains('Pre')`; the hybrid labels renamed Glut-Pre-1/2/3 to Glut-NSC/NB/Im (snrna/naming/hybrid_division_naming.qmd), so a
substring match would silently keep them. The list is explicit and asserted instead.

    adata.h5ad   .raw = raw integer counts (what SCENIC+ reads); X = log-normalized, HVG, scaled (display only)

    python make_glut_adata.py
"""
import os

import anndata as ad
import numpy as np
import scanpy as sc

B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/")
SRC = B + "ra-arco-hvc-nc_hybrid/data/adata.h5ad"          # X = counts, cluster = cluster_hybrid
OUT_DIR = B + "ra-arco-hvc-nc_glut_hybrid/anndata_rna/"
PRECURSORS = ["Glut-NSC", "Glut-NB", "Glut-Im"]

adata = sc.read_h5ad(SRC)
assert np.all(adata.X.data == np.round(adata.X.data)), "source X is not integer counts"

cl = adata.obs["cluster"].astype(str)
keep = cl.str.contains("Glut") & ~cl.isin(PRECURSORS)
adata = adata[keep.values].copy()
assert not adata.obs["cluster"].astype(str).str.contains("Pre").any()
assert set(PRECURSORS).isdisjoint(adata.obs["cluster"].astype(str)), "precursor cluster survived the filter"
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
