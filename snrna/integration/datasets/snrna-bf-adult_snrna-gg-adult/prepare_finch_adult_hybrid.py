"""Prepare the UPDATED adult finch (song-system-grn, hybrid-labeled) h5ad for SATURN/SAMap/GSI/CCA.

Source: group/sc_datasets/song-system-grn/obj_hybrid_labels.h5ad -- same 34,295 cells x
17,029 genes as the original bf_adult.h5ad (identical gene ID scheme, verified X already
raw integer counts), but with a refined 'celltype_hybrid' label column (47 distinct types)
replacing the flat 'cluster' scheme used in the original run. 659 cells have no
celltype_hybrid assignment and are dropped.

Notably, celltype_hybrid's naming already encodes the song/non-song structure directly:
position-restricted song-nucleus Glut populations get their own dedicated clusters
(Glut-DACH2-HVCra, -HVCra-Pre, -HVCx, Glut-CACNA1H-RA, Glut-DACH2-LMANco), while
Glut-DACH2-LMANsh (LMAN shell, non-song per this project's prior determination) and the
widely-shared non-song populations (Glut-DACH2-1..8, Glut-CACNA1H-1..4, spanning
NC/NR/Arco/HVC diffusely) are NOT further subdivided by anatomical position -- consistent
with the manuscript's claim that only the song nuclei contain position-specific derived
types. New developmental/ambiguous populations not in the original scheme (Glut-Im, -NB,
-NSC, -GABA, Glut-SATB2-1) are kept as-is; no attempt is made here to fold them into
song/non-song, since they are not clearly analogous to the adult excitatory types being
compared. That grouping decision is for the downstream composite/statistical analysis,
not this data-prep step.

Written to a NEW path (bf_adult_hybrid.h5ad) so as not to clobber the original
bf_adult.h5ad or any of its downstream results (train/saturn_results/macro2000_hv8000_seed0,
analysis/macro2000_hv8000_seed0, samap results, gsi_corr_gg_adult.csv, cca gg_adult_* --
all keyed to the OLD label scheme and left untouched).
"""
from pathlib import Path

import anndata as ad
import numpy as np

IN_H5AD = "/private/home/bcolquit/group/sc_datasets/song-system-grn/obj_hybrid_labels.h5ad"
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data/bf_adult_hybrid.h5ad")


def main():
    a = ad.read_h5ad(IN_H5AD)
    print(f"loaded {a.shape}")

    import scipy.sparse as sp
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    dev = np.abs(X.data - np.round(X.data)).max()
    print(f"max deviation from integer: {dev:.2e}")
    assert dev < 1e-3, "X is not already raw counts -- check source file"
    X = X.astype(np.float32)

    n_missing = a.obs["celltype_hybrid"].isna().sum()
    print(f"dropping {n_missing} cells with no celltype_hybrid assignment")
    keep = a.obs["celltype_hybrid"].notna().values

    out = ad.AnnData(X=X[keep], obs=a.obs.loc[keep].copy(), var=a.var.copy())
    out.obs["cluster"] = out.obs["celltype_hybrid"].astype(str)
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
