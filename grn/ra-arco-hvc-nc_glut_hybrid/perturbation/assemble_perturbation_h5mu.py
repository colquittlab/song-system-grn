#!/usr/bin/env python3
"""Perturbation-ready scplusmdata for a SCENIC+ config that never built the full ~35 GB merged file (build_scplus_mudata: false).

perturbation_parallel_template.ipynb needs the real expression matrix (`scRNA_counts`) next to the AUC matrices and eRegulon tables. In a full
scplusmdata.h5mu that modality is the GEX anndata's `.raw` (raw integer counts; `use_raw_for_GEX_anndata`) restricted to the cells present in
both assays, with the GEX obs / var / obsm (scenicplus.data_wrangling.process_multiome_data) -- checked against the all-cell hybrid config3's full file
(16,467 x 15,948, integer-valued). This rebuilds exactly that modality from anndata_rna/adata.h5ad, keeps SCENIC+'s own ScenicPlusMuData layout,
and replaces only the accessibility modality (never read here) with a 1-column placeholder, as assemble_slim_scplusmdata.py does.

    python assemble_perturbation_h5mu.py <config outs dir> [out.h5mu] [--lognorm]   # default: <outs>/scplusmdata_perturb.h5mu
With --lognorm the expression is log1p(CP10k) of those counts (scanpy normalize_total target_sum=1e4, then log1p) and the default name is
scplusmdata_perturb_lognorm.h5mu; this is NOT what SCENIC+ saw when it inferred the eRegulons (raw counts), so use it only to compare scales.
Inputs in <outs>: AUCell_direct.h5mu AUCell_extended.h5mu eRegulon_direct.tsv eRegulons_extended.tsv.  ADATA (env) overrides the expression source.
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

B = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
ADATA = os.environ.get("ADATA", B + "ra-arco-hvc-nc_glut_hybrid/anndata_rna/adata.h5ad")
LOGNORM = "--lognorm" in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith("--")]
outs = Path(args[0])
out = Path(args[1]) if len(args) > 1 else outs / ("scplusmdata_perturb_lognorm.h5mu" if LOGNORM else "scplusmdata_perturb.h5mu")

auc_d = mudata.read(str(outs / "AUCell_direct.h5mu"))
auc_e = mudata.read(str(outs / "AUCell_extended.h5mu"))
meta_d = pd.read_table(outs / "eRegulon_direct.tsv")
meta_e = pd.read_table(outs / "eRegulons_extended.tsv")
cells = list(auc_d["Gene_based"].obs_names)
assert cells == list(auc_e["Gene_based"].obs_names), "cell order differs between AUC files"

full = ad.read_h5ad(ADATA)
missing = [c for c in cells if c not in full.obs_names]
assert not missing, f"{len(missing)} cells absent from {ADATA}, e.g. {missing[:3]}"
raw = full.raw.to_adata()                                  # raw integer counts, all genes
raw = raw[cells].copy()
assert np.all(raw.X.data == np.round(raw.X.data)), "expression is not integer counts"
if LOGNORM:
    import scanpy as sc
    raw.X = csr_matrix(raw.X, dtype=np.float32)
    sc.pp.normalize_total(raw, target_sum=1e4)
    sc.pp.log1p(raw)
rna = ad.AnnData(X=csr_matrix(raw.X, dtype=np.float32), obs=full.obs.loc[cells].copy(), var=raw.var.copy(),
                 obsm={k: np.asarray(v)[[full.obs_names.get_loc(c) for c in cells]] for k, v in full.obsm.items()})
for c in rna.obs.columns:
    if rna.obs[c].dtype == object:
        rna.obs[c] = rna.obs[c].astype(str).astype("category")
atac = ad.AnnData(X=csr_matrix((len(cells), 1), dtype=np.float32), obs=rna.obs.copy(), var=pd.DataFrame(index=["slim_placeholder"]))

md = ScenicPlusMuData(acc_gex_mdata=mudata.MuData({"scRNA": rna, "scATAC": atac}),
                      e_regulon_auc_direct=auc_d, e_regulon_auc_extended=auc_e,
                      e_regulon_metadata_direct=meta_d, e_regulon_metadata_extended=meta_e)
md.write_h5mu(str(out))
print(f"wrote {out} ({out.stat().st_size / 1e6:.0f} MB): {len(cells)} cells x {rna.n_vars} genes ({'log1p CP10k' if LOGNORM else 'raw counts'}); "
      f"{auc_d['Gene_based'].n_vars} direct and {auc_e['Gene_based'].n_vars} extended eRegulons")
