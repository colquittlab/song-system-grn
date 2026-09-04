"""Prepare the adult finch (song-system-grn, Bengalese finch, lonStrDom2) h5ad for SATURN.

Unlike bf_dev.h5ad (log1p'd, needed expm1 recovery), this dataset's X is already raw
integer counts (verified: 0 deviation from integer on a 500k-entry sample). var_names
are lonStrDom2 (GCF_005870125.1) LOC-style RefSeq gene IDs, matched exactly against the
existing lonStrDom2 ESM1b embedding (85.4% coverage) -- no new embedding needed.

obs only carries 'cluster' (50 song-system cell types) and 'position' (6 song nuclei:
hvc, nr, nc, arco, lman, ra). All clusters hold >=74 cells, so no rare-label coalescing
is needed (contrast the La Manno mouse prep, which had a long single-digit tail).
"""
from pathlib import Path

import anndata as ad
import numpy as np

IN_H5AD = "/private/home/bcolquit/group/sc_datasets/song-system-grn/obj_clustered.h5ad/obj_clustered.h5ad"
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data/bf_adult.h5ad")


def main():
    a = ad.read_h5ad(IN_H5AD)
    print(f"loaded {a.shape}")

    import scipy.sparse as sp
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    dev = np.abs(X.data - np.round(X.data)).max()
    print(f"max deviation from integer: {dev:.2e}")
    assert dev < 1e-3, "X is not already raw counts -- check source file"
    X = X.astype(np.float32)

    out = ad.AnnData(X=X, obs=a.obs.copy(), var=a.var.copy())
    out.obs["cluster"] = out.obs["cluster"].astype(str)
    out.obs["position"] = out.obs["position"].astype(str)
    out.obs["species"] = "finch"
    out.var_names_make_unique()

    vc = out.obs["cluster"].value_counts()
    print(f"labels: {len(vc)}  median {int(vc.median())}  min {int(vc.min())}")
    print(f"median UMIs/cell: {np.median(np.asarray(out.X.sum(1)).ravel()):.0f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(OUT, compression="gzip")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
