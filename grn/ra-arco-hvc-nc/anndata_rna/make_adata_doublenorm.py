#!/usr/bin/env python3
"""Build the LEGACY, deliberately wrong SCENIC+ expression input, to test what the earlier run's mistake did.

The earlier run's `.raw` held log1p(normalize_total(SCT data)): the exported X was already SCT log data, and
`adata.raw = adata` shares the sparse matrix, so the in-place normalize_total/log1p that followed also rewrote
`.raw`. SCENIC+ reads `.raw` (prepare_GEX_ACC), so that matrix was the GBM input. This script recreates that
matrix explicitly (no reliance on the aliasing) on top of the production adata.h5ad, so cells, labels, barcodes
and gene order are identical and only `.raw` differs:

    adata.h5ad              .raw = raw integer counts            (correct; what SCENIC+ expects)
    adata_doublenorm.h5ad   .raw = log1p(normalize_total(SCT data))   (legacy; used by config12/13 only)

Never use adata_doublenorm.h5ad for real analyses.

    python make_adata_doublenorm.py
"""
import anndata as ad
import numpy as np
import scanpy as sc

B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/")
prod = ad.read_h5ad(B + "anndata_rna/adata.h5ad")        # raw = counts
sct = ad.read_h5ad(B + "data/adata_sct_data.h5ad")       # X = SCT log data, same cell order (export_h5ad_hybrid.R)

# same cells in the same order: labels, and barcodes once the notebook's "-N" -> "-1" rewrite is applied
assert prod.n_obs == sct.n_obs, (prod.n_obs, sct.n_obs)
for col in ("cluster", "position", "assignment", "sample_short", "run", "replicate"):
    assert (prod.obs[col].astype(str).values == sct.obs[col].astype(str).values).all(), col
bc_prod = prod.obs_names.str.split("___").str[0]
bc_sct = sct.obs_names.str.replace("-[0-9]", "-1", regex=True)
assert (bc_prod == bc_sct).all(), "barcode order differs"
assert list(prod.raw.var_names) == list(sct.var_names), "gene order differs"

# the legacy matrix: normalize SCT log data to 1e4 per cell, then log1p again
tmp = ad.AnnData(X=sct.X.copy().astype(np.float32), obs=prod.obs[[]].copy(), var=sct.var[[]].copy())
sc.pp.normalize_total(tmp, target_sum=1e4)
sc.pp.log1p(tmp)

out = prod.copy()
out.raw = ad.AnnData(X=tmp.X, obs=prod.obs.copy(), var=prod.raw.var.copy())
out.uns["raw_matrix"] = "LEGACY log1p(normalize_total(SCT data)); not raw counts. Test input only."

x = out.raw.X.data
print("raw (legacy):", out.raw.shape, "integer-valued:", bool(np.all(x == np.round(x))), "max: %.3f" % x.max(),
      "median per-cell sum: %.0f" % np.median(np.asarray(out.raw.X.sum(1)).ravel()))
c = prod.raw.X.data
print("raw (production counts):", prod.raw.shape, "integer-valued:", bool(np.all(c == np.round(c))), "max: %.0f" % c.max())
assert not np.all(x == np.round(x)) and np.all(c == np.round(c))

# Check against the earlier run's actual h5ad when it is reachable: on the cells both runs share, the rebuilt `.raw`
# must equal the old one up to float32 rounding (measured 1.5e-5 max abs diff over all 16,588 shared cells).
OLD = ("/mnt/nest/bcolquit/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
       "ra-arco-hvc-nc_seurat-snrna-clustering/anndata_rna/adata.h5ad")
try:
    old = ad.read_h5ad(OLD)
except OSError:
    old = None
    print("earlier run's adata.h5ad not reachable; reconstruction not cross-checked")
if old is not None:
    shared = out.obs_names.intersection(old.obs_names)
    assert list(old.raw.var_names) == list(out.raw.var_names), "gene order differs from the earlier run"
    diff = old.raw[shared].X - out.raw[shared].X
    worst = abs(diff.data).max()
    print(f"cross-check vs earlier run's .raw on {len(shared)} shared cells: max abs diff {worst:.1e}")
    assert len(shared) == out.n_obs and worst < 1e-3, "rebuilt legacy .raw does not match the earlier run's"

out.write_h5ad(B + "anndata_rna/adata_doublenorm.h5ad")
print("wrote", B + "anndata_rna/adata_doublenorm.h5ad")
