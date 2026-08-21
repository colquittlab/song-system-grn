"""Export the per-species inputs SCTransform needs for the SCT-normed GSI variant of the
glut-only chicken hybrid-label comparison.

Replicates run_gsi.py's prep() pipeline (robust_genes on the COMPLETE data -> subsample
cells to median population size -> thin UMIs in large populations) EXACTLY, using the
same zaremba_prep functions and the same ortholog file/seed, but stops BEFORE the
CPM library-size-normalization step -- SCTransform (run separately, per species, by
R/sct_transform_export.R in finch-integration-toolkit) takes over that role instead.

Orthologs: OrthoFinder (not RBH) for this specific finch x chicken pair, matching
assemble_glutonly.py's documented convention ("per explicit request... matches Zaremba's
own cited method with no substitution") -- NOT the RBH set used by the RPCA/SAMap side of
this analysis. GSI and RPCA have always used different ortholog sources in this project.
"""
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from zaremba_prep import robust_genes, subsample_cells, subsample_umis

GLUTONLY = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")
ORTHOLOGS = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/data/orthologs_orthofinder_bfgg_adult_aliased.csv")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/data/seurat_sct_gsi")
SEED = 0

QUERY_H5AD = GLUTONLY / "bf_adult_glut_hybrid.h5ad"
REF_H5AD = GLUTONLY / "gg_adult_ex.h5ad"
QUERY_LABEL = "cluster"
REF_LABEL = "cluster"


def export_species(h5ad: Path, label_col: str, genes: list, name: str, out_prefix: Path):
    print(f"\n--- {name} ---")
    a = ad.read_h5ad(h5ad)
    a.var_names = pd.Index([str(g) for g in a.var_names])
    a = a[:, genes].copy()
    a.obs[label_col] = a.obs[label_col].astype(str)
    labs_full = a.obs[label_col].values
    ok = robust_genes(a, label_col)
    keep = subsample_cells(labs_full, SEED)
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    X = X[keep]
    labs = labs_full[keep]
    X = subsample_umis(X, labs, SEED)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    mmwrite(str(out_prefix) + "_matrix.mtx", X.T.tocoo())  # genes x cells, matching export_for_seurat.py
    pd.Series(a.var_names).to_csv(str(out_prefix) + "_genes.tsv", index=False, header=False)
    cells = pd.Index(a.obs_names[keep])
    pd.Series(cells).to_csv(str(out_prefix) + "_cells.tsv", index=False, header=False)
    pd.DataFrame({"label": labs}, index=cells).to_csv(str(out_prefix) + "_labels.csv")
    print(f"  wrote {X.shape[1]:,} cells x {X.shape[0]:,} genes -> {out_prefix}_*")
    return ok


def main():
    orth = pd.read_csv(ORTHOLOGS)
    print(f"orthologs (1:1, OrthoFinder): {len(orth):,}")
    q = ad.read_h5ad(QUERY_H5AD, backed="r"); r = ad.read_h5ad(REF_H5AD, backed="r")
    orth = orth[orth.finch.isin(set(q.var_names)) & orth.mouse.isin(set(r.var_names))]
    print(f"  present in both datasets: {len(orth):,}")
    orth.to_csv(OUT / "orthologs_used.csv", index=False)

    okq = export_species(QUERY_H5AD, QUERY_LABEL, list(orth.finch), "finch (query)", OUT / "finch")
    okr = export_species(REF_H5AD, REF_LABEL, list(orth.mouse), "chicken (reference)", OUT / "chicken")

    both = okq & okr
    print(f"\ngenes robust in BOTH species: {both.sum():,}/{len(both):,}")
    np.save(OUT / "both_robust_mask.npy", both)
    print(f"wrote {OUT}/both_robust_mask.npy (apply AFTER SCTransform, in the same ortholog row order)")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()
