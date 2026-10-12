#!/usr/bin/env python3
"""Pseudobulk accessibility (fragments per million, CPM) of ALL consensus regions in the four cell types behind the two contrasts: RA and C1H-1
(RA vs C1H-1), PVALB-1 and LAMP5 (PVALB-1 vs LAMP5). Source: the fragment matrix of the 40-topic cisTopic object (499,348 regions, the consensus of
configs 34-39). CPM = fragments in the region / total fragments in all regions of that cell type x 1e6, so values are comparable across cell types of
different size and depth. Read by scenicplus_regulon_region_heatmaps.R, which plots linear differences (not ratios) and the accessibility itself.

    python scenicplus_region_accessibility_all_regions.py
"""
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

STORE = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid")
OUT = Path("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/config37/regulon_region_heatmaps").expanduser()
OUT.mkdir(parents=True, exist_ok=True)
GROUPS = {"RA": "Glut-CACNA1H-RA", "C1H1": "Glut-CACNA1H-1", "PV1": "GABA-MGE-PVALB-1", "LAMP5": "GABA-MGE-LAMP5"}

o = pickle.load(open(STORE / "pycisTopic/cistopic_obj_glut_k40.pkl", "rb"))
cl = o.cell_data.loc[o.cell_names, "cluster"].astype(str).values
print({k: int((cl == v).sum()) for k, v in GROUPS.items()}, "cells")
fm = o.fragment_matrix.tocsr()
idx = {v: i for i, v in enumerate(GROUPS.values())}
onehot = sp.csr_matrix((np.ones(len(cl)), (np.arange(len(cl)), [idx.get(c, len(idx)) for c in cl])), shape=(len(cl), len(idx) + 1))[:, : len(idx)]
counts = (fm @ onehot).toarray().astype(np.float64)                 # regions x 4
cpm = counts / counts.sum(0, keepdims=True) * 1e6
df = pd.DataFrame(cpm, index=o.region_names, columns=list(GROUPS)).round(4)
df.index.name = "region"
df.to_csv(OUT / "region_accessibility_cpm.tsv.gz", sep="\t")
print("regions:", len(df), "| median CPM per group:", df.median().round(3).to_dict(), "| 99th percentile:", df.quantile(0.99).round(1).to_dict())
