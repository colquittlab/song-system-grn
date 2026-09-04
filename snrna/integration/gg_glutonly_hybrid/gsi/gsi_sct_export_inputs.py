"""Export the per-species inputs SCTransform needs for the SCT-normed GSI variant of the
glut-only chicken hybrid-label analysis -- generalized to run across all four finch/chicken
subset combinations now in play in this analysis, so GSI stays in sync with the same data
subsets used by the RPCA/SAMap side:

  1. base:                finch=bf_adult_glut_hybrid.h5ad,          chicken=gg_adult_ex.h5ad
  2. plusSATB2:            finch=bf_adult_glut_hybrid_plusSATB2.h5ad, chicken=gg_adult_ex.h5ad
  3. noMeso:               finch=bf_adult_glut_hybrid.h5ad,          chicken=gg_adult_ex_noMeso.h5ad
  4. plusSATB2_noPre:      finch=bf_adult_glut_hybrid_plusSATB2.h5ad, chicken=gg_adult_ex_noPre.h5ad

Replicates run_gsi.py's prep() pipeline (robust_genes on the COMPLETE data -> subsample
cells to median population size -> thin UMIs in large populations) EXACTLY per subset,
using the same zaremba_prep functions and the same ortholog file/seed, but stops BEFORE
the CPM library-size-normalization step -- SCTransform (run separately, per species, by
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
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/gsi/data")
SEED = 0
QUERY_LABEL = "cluster"
REF_LABEL = "cluster"

TESTS = [
    ("base", "bf_adult_glut_hybrid.h5ad", "gg_adult_ex.h5ad"),
    ("plusSATB2", "bf_adult_glut_hybrid_plusSATB2.h5ad", "gg_adult_ex.h5ad"),
    ("noMeso", "bf_adult_glut_hybrid.h5ad", "gg_adult_ex_noMeso.h5ad"),
    ("plusSATB2_noPre", "bf_adult_glut_hybrid_plusSATB2.h5ad", "gg_adult_ex_noPre.h5ad"),
]


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
    OUT.mkdir(parents=True, exist_ok=True)
    orth_all = pd.read_csv(ORTHOLOGS)
    print(f"orthologs (1:1, OrthoFinder): {len(orth_all):,}")

    for tag, query_name, ref_name in TESTS:
        print(f"\n=== {tag} ===")
        query_h5ad = GLUTONLY / query_name
        ref_h5ad = GLUTONLY / ref_name
        q = ad.read_h5ad(query_h5ad, backed="r")
        r = ad.read_h5ad(ref_h5ad, backed="r")
        orth = orth_all[orth_all.finch.isin(set(q.var_names)) & orth_all.mouse.isin(set(r.var_names))]
        print(f"  present in both datasets: {len(orth):,}")
        orth.to_csv(OUT / f"orthologs_used_{tag}.csv", index=False)

        okq = export_species(query_h5ad, QUERY_LABEL, list(orth.finch), f"finch ({tag})", OUT / f"finch_{tag}")
        okr = export_species(ref_h5ad, REF_LABEL, list(orth.mouse), f"chicken ({tag})", OUT / f"chicken_{tag}")

        both = okq & okr
        print(f"  genes robust in BOTH species: {both.sum():,}/{len(both):,}")
        np.save(OUT / f"both_robust_mask_{tag}.npy", both)


if __name__ == "__main__":
    main()
