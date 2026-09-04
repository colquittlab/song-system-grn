"""Reassemble SCT-corrected h5ads for SAMap from R/sct_transform_export.R's mtx output +
the original per-species .obs metadata (cell order/identity unaffected by SCTransform;
gene set may have shrunk to SCTransform's own filter, which is fine -- SAMap works with
whatever genes are present in each object's X)."""
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io as sio

IN = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/data")
GLUTONLY = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")

SPECIES = [("finch", GLUTONLY / "bf_adult_glut_hybrid.h5ad"), ("chicken", GLUTONLY / "gg_adult_ex.h5ad")]


def assemble(name: str, orig_h5ad: Path):
    prefix = IN / f"{name}_sct"
    X = sio.mmread(str(prefix) + "_matrix.mtx").tocsr().T.tocsr()  # -> cells x genes
    genes = pd.read_csv(str(prefix) + "_genes.tsv", header=None)[0].values
    cells = pd.read_csv(str(prefix) + "_cells.tsv", header=None)[0].values

    orig = ad.read_h5ad(orig_h5ad, backed="r")
    obs = orig.obs.loc[cells].copy()
    assert list(cells) == list(obs.index), "cell order mismatch"

    out = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    out.var_names_make_unique()
    out_path = IN / f"{name}_sct.h5ad"
    out.write_h5ad(out_path, compression="gzip")
    print(f"{name}: {out.shape} -> {out_path}")


if __name__ == "__main__":
    for name, orig in SPECIES:
        assemble(name, orig)
