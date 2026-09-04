"""Label-permutation significance test for the GSI-SCT finch-chicken correlation
matrices, across all four finch/chicken subsets (base, plusSATB2, noMeso,
plusSATB2_noPre).

Null model: for each finch cluster, is its observed max GSI correlation (across all
chicken clusters) higher than what a RANDOM group of the same number of finch cells would
achieve? Concretely: shuffle the per-cell finch cluster-label vector (reassigning which
cells get which label, preserving each cluster's exact size), recompute the finch
GSI-Z profile from scratch under the shuffled labels, correlate against the REAL
(unpermuted) chicken GSI-Z profile, and record each label's max-r under that shuffle.
Repeat N_PERM times to build an empirical null distribution per finch cluster (pooling
over all clusters of that size across permutations would be an alternative but is not
needed here -- shuffling the actual label vector already gives each real cluster its own
same-size null draw every permutation, since label identity and size stay linked).
Chicken labels are NOT permuted -- the question is whether THIS finch cluster's own
cells support its specificity profile being real, not whether the chicken side is
structured at all (already well established).

p-value = (1 + #{null max-r >= observed max-r}) / (N_PERM + 1), one-sided (is the
observed match stronger than chance), with a Benjamini-Hochberg FDR correction across
the ~18-19 finch clusters tested per subset.

Mirrors run_gsi_sct.py's exact GSI/rank/z-score computation (zaremba_prep.gsi, then
per-species rank + z-score, then Zq.T @ Zr / n_shared_genes) so permutation p-values are
directly comparable to that script's saved correlation matrices -- including replicating
its shared = min(len(Gq), len(Gr)) gene-count truncation, recomputed fresh each
permutation since gsi()'s own nonzero-total gene filter can vary slightly with the
relabeling.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy.io as sio

sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from zaremba_prep import gsi

IN = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/gsi/data")
OUT_DIR = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/gsi/results")
N_PERM = 2000
SEED = 0

TESTS = ["base", "plusSATB2", "noMeso", "plusSATB2_noPre"]


def load_percell(prefix: Path, ortholog_genes: list):
    """Per-cell SCT-corrected counts, reindexed to the ortholog gene order (missing
    genes -> 0 row, matching run_gsi_sct.py's load_sct_means)."""
    X = sio.mmread(str(prefix) + "_matrix.mtx").tocsr()  # genes x cells, SCT's own gene set/order
    genes = pd.read_csv(str(prefix) + "_genes.tsv", header=None)[0].values
    cells = pd.read_csv(str(prefix) + "_cells.tsv", header=None)[0].values
    gene_pos = {g: i for i, g in enumerate(genes)}
    idx = np.array([gene_pos.get(g, -1) for g in ortholog_genes])
    out = np.zeros((len(ortholog_genes), X.shape[1]), dtype=np.float64)
    ok = idx >= 0
    out[ok, :] = X[idx[ok], :].toarray()
    return out, cells


def means_by_label(Xg: np.ndarray, labels: np.ndarray, cats):
    """One-hot matmul instead of a per-cluster boolean-mask loop -- ~50x faster
    (BLAS-optimized single matrix multiply vs. len(cats) separate column-slice + mean
    operations), which matters here since this runs inside the permutation loop."""
    cat_idx = {c: i for i, c in enumerate(cats)}
    codes = np.fromiter((cat_idx[l] for l in labels), dtype=np.int64, count=len(labels))
    I = np.zeros((len(labels), len(cats)))
    I[np.arange(len(labels)), codes] = 1.0
    sums = Xg @ I
    counts = I.sum(0)
    return pd.DataFrame(sums / counts, columns=cats)


def to_Z(mean_df: pd.DataFrame):
    """gsi() -> rank -> z-score per zaremba_prep/run_gsi_sct.py's convention. gsi() only
    drops ROWS (genes with zero total), never columns, so the cluster/column set and
    order are always preserved."""
    G = gsi(mean_df)
    Rk = G.rank(axis=0).values
    Z = (Rk - Rk.mean(0)) / np.where(Rk.std(0) == 0, 1, Rk.std(0))
    return Z


def run_one(tag: str) -> pd.DataFrame:
    print(f"\n=== {tag} ===")
    orth = pd.read_csv(IN / f"orthologs_used_{tag}.csv")
    both = np.load(IN / f"both_robust_mask_{tag}.npy")

    Xq_full, q_cells = load_percell(IN / f"finch_{tag}_sct", list(orth.finch))
    Xr_full, r_cells = load_percell(IN / f"chicken_{tag}_sct", list(orth.mouse))
    Xq = Xq_full[both]
    Xr = Xr_full[both]

    q_labs = pd.read_csv(IN / f"finch_{tag}_labels.csv", index_col=0)["label"].reindex(q_cells).values
    r_labs = pd.read_csv(IN / f"chicken_{tag}_labels.csv", index_col=0)["label"].reindex(r_cells).values
    assert not pd.isna(q_labs).any() and not pd.isna(r_labs).any(), "cell order mismatch"
    q_cats = pd.unique(q_labs)
    r_cats = pd.unique(r_labs)
    print(f"finch: {Xq.shape[1]:,} cells, {len(q_cats)} clusters | chicken: {Xr.shape[1]:,} cells, {len(r_cats)} clusters")

    Mr = means_by_label(Xr, r_labs, r_cats)
    Zr = to_Z(Mr)

    Mq_real = means_by_label(Xq, q_labs, q_cats)
    Zq_real = to_Z(Mq_real)
    shared = min(Zq_real.shape[0], Zr.shape[0])
    C_real = Zq_real[:shared].T @ Zr[:shared] / shared
    real_max = pd.Series(C_real.max(axis=1), index=q_cats)

    rng = np.random.default_rng(SEED)
    null_max = {c: np.empty(N_PERM) for c in q_cats}
    for p in range(N_PERM):
        perm_labels = rng.permutation(q_labs)
        Mq_p = means_by_label(Xq, perm_labels, q_cats)
        Zq_p = to_Z(Mq_p)
        shared_p = min(Zq_p.shape[0], Zr.shape[0])
        Cp = Zq_p[:shared_p].T @ Zr[:shared_p] / shared_p
        row_max = Cp.max(axis=1)
        for i, c in enumerate(q_cats):
            null_max[c][p] = row_max[i]
        if (p + 1) % 500 == 0:
            print(f"  {p + 1}/{N_PERM} permutations")

    rows = []
    for c in q_cats:
        arr = null_max[c]
        obs = real_max[c]
        p_val = (1 + (arr >= obs).sum()) / (1 + len(arr))
        rows.append({
            "finch_cluster": c,
            "observed_max_r": round(float(obs), 4),
            "null_mean": round(float(arr.mean()), 4),
            "null_p95": round(float(np.percentile(arr, 95)), 4),
            "p_value": round(float(p_val), 4),
        })
    res = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    m = len(res)
    bh = (res["p_value"] * m / (res.index + 1)).clip(upper=1.0)
    res["fdr_bh"] = bh[::-1].cummin()[::-1].round(4)

    out_csv = OUT_DIR / f"gsi_sct_permutation_{tag}.csv"
    res.to_csv(out_csv, index=False)
    pd.set_option("display.width", 160)
    print(res.to_string(index=False))
    print(f"wrote {out_csv}")
    return res


if __name__ == "__main__":
    for tag in TESTS:
        run_one(tag)
