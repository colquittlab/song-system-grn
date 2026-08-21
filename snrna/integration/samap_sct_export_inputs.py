"""Export the FULL per-species raw-count matrices (all cells, all genes -- no ortholog
restriction, no subsampling) for SCTransform, as input to the SAMap-on-SCT variant of the
glut-only chicken hybrid-label comparison.

Unlike the GSI-SCT export, SAMap does NOT restrict to a shared ortholog gene space before
normalization -- it uses its own BLAST-graph homology across each species' FULL
transcriptome, and never subsamples cells (subsampling is a GSI/label-transfer-specific
step in this project's methods, not part of SAMap's own algorithm). So SCTransform here
runs on the complete per-species objects already used by run_samap.py.
"""
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

GLUTONLY = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/data/seurat_sct_samap")

SPECIES = [
    ("finch", GLUTONLY / "bf_adult_glut_hybrid.h5ad"),
    ("chicken", GLUTONLY / "gg_adult_ex.h5ad"),
]


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
    OUT.mkdir(parents=True, exist_ok=True)
    for name, h5ad in SPECIES:
        export(name, h5ad)
