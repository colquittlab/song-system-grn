"""Export inputs for a MetaNeighborUS pilot on the glut-only finch x chicken excitatory
comparison (plusSATB2 finch subset x noPre chicken reference -- the current/latest test
set, so this pilot is directly comparable to the RPCA/SAMap results already produced for
it).

Rationale: unlike Seurat's TransferData (a nearest-neighbour VOTE, always assigns some
predicted.id even when no real match exists -- the artifact reciprocal_score.py's
sqrt(fwd*rev) symmetrization works around) and SAMap (reciprocal by construction, but
still a single best-hit call), MetaNeighborUS reports a continuous AUROC for EVERY
cluster pair. A cluster with no genuine cross-species counterpart should show uniformly
low (~0.5, chance-level) AUROC against everything, rather than being forced into a
"best available" match -- directly relevant to the HVCra question this whole analysis has
been chasing. This is a pilot: validate on a known-good non-song pair before trusting it
on HVCra, since MetaNeighbor is normally used within-species (its correlation-based
neighbour-voting network has not been validated at cross-species batch-effect scale here).

Same ortholog set (RBH, matching the RPCA/SAMap side of this analysis -- NOT the
OrthoFinder set used for GSI) and same population-balancing convention (subsample every
cluster down to the per-species MEDIAN population size, Zaremba-style) as the rest of
this pipeline -- but NO UMI thinning (that's GSI-specific) and no ortholog-set-specific
robust_genes filter beyond intersecting with what's present in both h5ads, since
MetaNeighborUS does its own highly-variable-gene selection internally (variableGenes()).

Output: a single combined (finch+chicken) log-normalized matrix on the shared ortholog
gene space, written as mtx + genes.tsv + cells.tsv + a meta.csv with 'species' (study_id)
and 'cluster' (cell_type) columns, for run_metaneighbor.R to assemble directly.
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite

import sys
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from zaremba_prep import subsample_cells

GLUTONLY = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")
ORTHOLOGS = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/data/orthologs_rbh_bfgg_adult_aliased.csv")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/metaneighbor/data")
SEED = 0

QUERY_H5AD = GLUTONLY / "bf_adult_glut_hybrid_plusSATB2.h5ad"
REF_H5AD = GLUTONLY / "gg_adult_ex_noPre.h5ad"


def main():
    orth = pd.read_csv(ORTHOLOGS)
    q = ad.read_h5ad(QUERY_H5AD, backed="r")
    r = ad.read_h5ad(REF_H5AD, backed="r")
    orth = orth[orth.finch.isin(set(q.var_names)) & orth.mouse.isin(set(r.var_names))]
    print(f"1:1 RBH orthologs present in both: {len(orth):,}")

    finch = ad.read_h5ad(QUERY_H5AD)[:, list(orth.finch)].copy()
    finch.var_names = pd.Index(list(orth.finch))
    chick = ad.read_h5ad(REF_H5AD)[:, list(orth.mouse)].copy()
    chick.var_names = pd.Index(list(orth.finch))  # rename to shared (finch) symbol space

    finch.obs["cluster"] = finch.obs["cluster"].astype(str)
    chick.obs["cluster"] = chick.obs["cluster"].astype(str)

    keep_f = subsample_cells(finch.obs["cluster"].values, SEED)
    keep_c = subsample_cells(chick.obs["cluster"].values, SEED)
    finch = finch[keep_f].copy()
    chick = chick[keep_c].copy()
    print(f"finch: {keep_f.sum():,} cells kept, {finch.obs['cluster'].nunique()} clusters")
    print(f"chicken: {keep_c.sum():,} cells kept, {chick.obs['cluster'].nunique()} clusters")

    finch.obs_names = [f"bf_{c}" for c in finch.obs_names]
    chick.obs_names = [f"gg_{c}" for c in chick.obs_names]
    finch.obs["species"] = "finch"
    chick.obs["species"] = "chicken"

    combo = ad.concat([finch, chick], join="outer")
    sc.pp.normalize_total(combo, target_sum=1e4)
    sc.pp.log1p(combo)

    OUT.mkdir(parents=True, exist_ok=True)
    X = combo.X.tocsr() if sp.issparse(combo.X) else sp.csr_matrix(combo.X)
    mmwrite(str(OUT / "combined_matrix.mtx"), X.T.tocoo())  # genes x cells
    pd.Series(combo.var_names).to_csv(OUT / "combined_genes.tsv", index=False, header=False)
    pd.Series(combo.obs_names).to_csv(OUT / "combined_cells.tsv", index=False, header=False)
    combo.obs[["species", "cluster"]].to_csv(OUT / "combined_meta.csv")
    print(f"\nwrote {X.shape[1]:,} genes x {X.shape[0]:,} cells -> {OUT}/combined_*")


if __name__ == "__main__":
    main()
