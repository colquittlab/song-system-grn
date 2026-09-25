"""Export the FULL per-species raw-count matrices (all cells, all genes -- no ortholog
restriction, no subsampling) for SCTransform, as input to the SAMap-on-SCT variant of the
glut-only chicken hybrid-label comparison.

Unlike the GSI-SCT export, SAMap does NOT restrict to a shared ortholog gene space before
normalization -- it uses its own BLAST-graph homology across each species' FULL
transcriptome, and never subsamples cells (subsampling is a GSI/label-transfer-specific
step in this project's methods, not part of SAMap's own algorithm). So SCTransform here
runs on the complete per-species objects already used by run_samap.py.

Generalized over the finch/chicken subset combinations used elsewhere in this analysis
(same table as gsi/gsi_sct_export_inputs.py); `--subset base` reproduces the original
finch_/chicken_ export, `--subset plusSATB2_noPre` writes finch_plusSATB2_/chicken_noPre_.
SCTransform is then run per species on EXACTLY that cell set (its regularized NB model is
fit on the cells present, so a subset is re-normalized rather than sliced from the base
SCT output).
"""
import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

GLUTONLY = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/data")

SUBSETS = {   # subset -> [(export name, source h5ad), ...]
    "base": [("finch", "bf_adult_glut_hybrid.h5ad"), ("chicken", "gg_adult_ex.h5ad")],
    "plusSATB2": [("finch_plusSATB2", "bf_adult_glut_hybrid_plusSATB2.h5ad"), ("chicken", "gg_adult_ex.h5ad")],
    "noMeso": [("finch", "bf_adult_glut_hybrid.h5ad"), ("chicken_noMeso", "gg_adult_ex_noMeso.h5ad")],
    "plusSATB2_noPre": [("finch_plusSATB2", "bf_adult_glut_hybrid_plusSATB2.h5ad"), ("chicken_noPre", "gg_adult_ex_noPre.h5ad")],
}


def export(name: str, h5ad: Path):
    print(f"\n--- {name} ---")
    a = ad.read_h5ad(h5ad)
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    out_prefix = OUT / name
    mmwrite(str(out_prefix) + "_matrix.mtx", X.T.tocoo())  # genes x cells
    pd.Series(a.var_names).to_csv(str(out_prefix) + "_genes.tsv", index=False, header=False)
    pd.Series(a.obs_names).to_csv(str(out_prefix) + "_cells.tsv", index=False, header=False)
    a.obs.to_csv(str(out_prefix) + "_obs.csv")
    print(f"  wrote {X.shape[0]:,} cells x {X.shape[1]:,} genes -> {out_prefix}_*")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", choices=list(SUBSETS), default="base")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for name, h5ad in SUBSETS[args.subset]:
        export(name, GLUTONLY / h5ad)
