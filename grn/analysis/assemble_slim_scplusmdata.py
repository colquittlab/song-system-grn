#!/usr/bin/env python3
"""Build scplusmdata_slim.h5mu for a SCENIC+ config from its small outputs, without the ~35 GB accessibility object.

`scenicplus grn_inference create_scplus_mudata` wraps the full accessibility+expression MuData together with the
AUCell results and eRegulon tables into scplusmdata.h5mu (~35 GB). Everything all.qmd reads from it is small: the four
AUC matrices, the eRegulon metadata, and the cells' cluster labels. This does the same wrapping with SCENIC+'s own
ScenicPlusMuData, but with 1-column placeholder scRNA/scATAC modalities that carry only the cell metadata, so the
layout is identical and the file is a few hundred MB. all.qmd falls back to the slim file when the full one is absent.

Inputs, all small and all produced by the workflow when `build_scplus_mudata: false`:
    <outs>/AUCell_direct.h5mu  AUCell_extended.h5mu  eRegulon_direct.tsv  eRegulons_extended.tsv
plus the production adata.h5ad for the cluster labels.

    python assemble_slim_scplusmdata.py <config outs dir> [out.h5mu]
"""
import os
import sys
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scenicplus.scenicplus_mudata import ScenicPlusMuData

B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/")
# Cluster labels come from the run's own adata.h5ad: ADATA env var, else the all-cell hybrid run; the glut run sets
# ADATA=.../ra-arco-hvc-nc_glut_hybrid/anndata_rna/adata.h5ad (see assemble_slim_all.sh, which takes ROOT).
ADATA = os.environ.get("ADATA", B + "ra-arco-hvc-nc_hybrid/anndata_rna/adata.h5ad")
OBS_COLS = ["assignment", "barcode", "cluster", "position", "replicate", "run", "sample_short"]  # as in the full file

outs = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else outs / "scplusmdata_slim.h5mu"

auc_d = mudata.read(str(outs / "AUCell_direct.h5mu"))
auc_e = mudata.read(str(outs / "AUCell_extended.h5mu"))
meta_d = pd.read_table(outs / "eRegulon_direct.tsv")
meta_e = pd.read_table(outs / "eRegulons_extended.tsv")

cells = list(auc_d["Gene_based"].obs_names)
assert cells == list(auc_e["Gene_based"].obs_names) == list(auc_d["Region_based"].obs_names), "cell order differs between AUC files"
full_obs = ad.read_h5ad(ADATA, backed="r").obs
missing = [c for c in cells if c not in full_obs.index]
assert not missing, f"{len(missing)} cells absent from adata.h5ad, e.g. {missing[:3]}"
obs = full_obs.loc[cells, OBS_COLS].copy()
for c in OBS_COLS:
    obs[c] = obs[c].astype(str).astype("category")  # categorical, as in the full file


def placeholder():
    return ad.AnnData(X=csr_matrix((len(cells), 1), dtype=np.float32), obs=obs.copy(), var=pd.DataFrame(index=["slim_placeholder"]))


slim = ScenicPlusMuData(
    acc_gex_mdata=mudata.MuData({"scRNA": placeholder(), "scATAC": placeholder()}),
    e_regulon_auc_direct=auc_d, e_regulon_auc_extended=auc_e,
    e_regulon_metadata_direct=meta_d, e_regulon_metadata_extended=meta_e)
slim.write_h5mu(str(out))
print(f"wrote {out} ({out.stat().st_size / 1e6:.0f} MB): {len(cells)} cells, "
      f"{auc_d['Gene_based'].n_vars} direct and {auc_e['Gene_based'].n_vars} extended eRegulons")
