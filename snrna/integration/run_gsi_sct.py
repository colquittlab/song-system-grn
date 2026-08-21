"""GSI correlation for the glut-only chicken hybrid-label comparison, SCT-normed variant
-- same structure as run_gsi.py's tail (mean-by-label -> GSI -> Spearman) but the
population means are computed from SCTransform's corrected counts (produced by
gsi_sct_export_inputs.py + R/sct_transform_export.R, run independently per species)
instead of run_gsi.py's CPM library-size-normalized counts. Cell/UMI subsampling upstream
of normalization is identical (same seed, same zaremba_prep functions) -- only the
normalization step differs, so this is a controlled comparison against the log-norm GSI
result, not a different pipeline end to end.

SCTransform applies its own gene filter, so a species' corrected-count matrix may be
missing some of the ortholog genes present in the pre-SCT export -- those are reindexed
in as all-zero rows (excluded again in practice by the robustness mask, which already
required nonzero expression in the COMPLETE pre-subsample data).
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy.io as sio

sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from zaremba_prep import gsi

IN = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/data/seurat_sct_gsi")
OUT_CSV = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/results/gsi_corr_gg_glutonly_hybrid_SCT.csv")


def load_sct_means(prefix: Path, ortholog_genes: list, labels_csv: Path) -> pd.DataFrame:
    """Load a species' SCT-corrected counts, reindex to the ortholog gene order (missing
    genes -> 0, since SCTransform's own filter can drop some), then mean-by-label."""
    X = sio.mmread(str(prefix) + "_matrix.mtx").tocsr()  # genes x cells, SCT's own gene set/order
    genes = pd.read_csv(str(prefix) + "_genes.tsv", header=None)[0].values
    cells = pd.read_csv(str(prefix) + "_cells.tsv", header=None)[0].values
    labs = pd.read_csv(labels_csv, index_col=0)["label"]
    labs = labs.reindex(cells).values
    assert not pd.isna(labs).any(), "cell order mismatch between SCT output and labels.csv"

    gene_pos = {g: i for i, g in enumerate(genes)}
    n_missing = sum(1 for g in ortholog_genes if g not in gene_pos)
    print(f"  {n_missing}/{len(ortholog_genes)} ortholog genes dropped by SCTransform's own filter -> 0")
    idx = [gene_pos.get(g, -1) for g in ortholog_genes]

    cats = pd.unique(labs)
    out = np.zeros((len(ortholog_genes), len(cats)), dtype=np.float64)
    for k, c in enumerate(cats):
        cols = labs == c
        sub = X[:, cols]
        means = np.asarray(sub.mean(1)).ravel()
        out[:, k] = [means[i] if i >= 0 else 0.0 for i in idx]
    return pd.DataFrame(out, columns=cats, index=ortholog_genes)


def main():
    orth = pd.read_csv(IN / "orthologs_used.csv")
    both = np.load(IN / "both_robust_mask.npy")
    print(f"genes robust in both species (from log-norm export step): {both.sum():,}/{len(both):,}")

    print("\n--- finch (query) ---")
    Mq = load_sct_means(IN / "finch_sct", list(orth.finch), IN / "finch_labels.csv")
    print("\n--- chicken (reference) ---")
    Mr = load_sct_means(IN / "chicken_sct", list(orth.mouse), IN / "chicken_labels.csv")

    Mq = Mq.loc[np.asarray(orth.finch)[both]]
    Mr = Mr.loc[np.asarray(orth.mouse)[both]]

    Gq, Gr = gsi(Mq), gsi(Mr)
    shared = min(len(Gq), len(Gr))
    Gq, Gr = Gq.iloc[:shared], Gr.iloc[:shared]

    Rq = Gq.rank(axis=0).values
    Rr = Gr.rank(axis=0).values
    Zq = (Rq - Rq.mean(0)) / np.where(Rq.std(0) == 0, 1, Rq.std(0))
    Zr = (Rr - Rr.mean(0)) / np.where(Rr.std(0) == 0, 1, Rr.std(0))
    C = pd.DataFrame(Zq.T @ Zr / Zq.shape[0], index=Gq.columns, columns=Gr.columns)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    C.round(5).to_csv(OUT_CSV)
    pd.set_option("display.width", 240)
    print(f"\nGSI (SCT-normed) Spearman matrix {C.shape} -> {OUT_CSV}")
    print(f"  range {C.values.min():+.3f} to {C.values.max():+.3f}, mean {C.values.mean():+.3f}")
    top = pd.DataFrame({"gsi_top": C.idxmax(axis=1), "r": C.max(axis=1).round(3)})
    print(top.to_string())


if __name__ == "__main__":
    main()
