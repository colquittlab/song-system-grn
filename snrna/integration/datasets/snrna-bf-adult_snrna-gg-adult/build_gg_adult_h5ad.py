"""Assemble the chicken developmental h5ad from the R/hdf5r export.

export_seurat.R writes the dgCMatrix slots directly (genes x cells, CSC) plus the
Seurat meta.data as CSV. This reconstructs a cells x genes AnnData with raw counts
in X, ready for SATURN's seurat_v3 HVG step.
"""
import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult")
DATA = BASE / "data"


def _dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


GENE_MAP = "/private/groups/colquittlab/saturn/embeddings/ens87_gg_geneid_to_symbol.tsv"


def main(label_col: list[str], min_label_cells: int, drop_labels: list[str]) -> None:
    with h5py.File(DATA / "gg_adult_counts.h5", "r") as f:
        shape = tuple(int(x) for x in f["shape"][:])          # (n_genes, n_cells)
        X = sp.csc_matrix(
            (f["data"][:].astype(np.float32), f["indices"][:], f["indptr"][:]),
            shape=shape,
        )
        genes = _dec(f["genes"][:])
        cells = _dec(f["cells"][:])

    print(f"loaded {shape[0]} genes x {shape[1]} cells, nnz={X.nnz:,}")
    X = X.T.tocsr()   # -> cells x genes

    md = pd.read_csv(DATA / "gg_adult_metadata.csv", low_memory=False)
    md = md.set_index("CellID").reindex(cells)
    print(f"metadata: {md.shape[0]} rows, {md.shape[1]} columns")

    a = ad.AnnData(X=X, obs=md, var=pd.DataFrame(index=pd.Index(genes, name=None)))

    # The matrix is indexed by Ensembl 87 gene IDs (ENSGALG00000*), plus ~28k XLOC
    # novel-transcript IDs from a StringTie annotation. The ESM1b embedding is keyed
    # by symbol against GRCg7b, whose IDs are a different series (ENSGALG00010*), so
    # translate through the Ensembl 87 GTF -- the same annotation behind EnsDb AH53208.
    sym = dict(l.split("\t") for l in open(GENE_MAP).read().strip().split("\n"))
    keep_g = [g in sym for g in a.var_names]
    print(f"genes: {a.n_vars} total, {sum(keep_g)} map to an Ensembl 87 symbol "
          f"({sum(g.startswith('XLOC') for g in a.var_names)} are XLOC novel transcripts)")
    a = a[:, keep_g].copy()
    a.var["ensembl_id_ens87"] = a.var_names
    a.var_names = pd.Index([sym[g] for g in a.var["ensembl_id_ens87"]])
    a.var_names_make_unique()

    # ~94.5% of values are already integral; the rest carry small fractional parts
    # (ambient correction / multimapped-read apportioning). Round so SATURN's
    # seurat_v3 HVG step receives count data.
    frac = float(np.mean(np.abs(a.X.data - np.round(a.X.data)) > 1e-6))
    print(f"non-integral entries: {frac:.1%} (max deviation "
          f"{np.abs(a.X.data - np.round(a.X.data)).max():.3f}) -- rounding")
    a.X.data = np.round(a.X.data).astype(np.float32)
    a.X.eliminate_zeros()

    lab = a.obs[label_col[0]].astype(str)
    for extra in label_col[1:]:
        lab = lab + " | " + a.obs[extra].astype(str)
    lab = lab.mask(lab.isin(["nan", "NA", "None", ""]), "Unassigned")
    if drop_labels:
        keep = ~lab.isin(drop_labels)
        print(f"dropping {(~keep).sum()} cells with labels {drop_labels}")
        a, lab = a[keep].copy(), lab[keep]

    counts = lab.value_counts()
    rare = counts[counts < min_label_cells].index
    if len(rare):
        print(f"merging {len(rare)} labels holding {int(counts[rare].sum())} cells -> 'Other (rare)'")
        lab = lab.mask(lab.isin(rare), "Other (rare)")

    a.obs["cluster"] = lab.values
    a.obs["species"] = "chicken"
    for extra_col in ["dissection", "anno_level_1", "class"]:
        if extra_col in a.obs.columns:
            a.obs[extra_col] = a.obs[extra_col].astype(str)

    vc = a.obs["cluster"].value_counts()
    print(f"final: {a.shape}  labels={len(vc)}  median={int(vc.median())}  min={int(vc.min())}")
    print(f"median UMIs/cell: {np.median(np.asarray(a.X.sum(1)).ravel()):.0f}")
    print(vc.head(15).to_string())

    out = DATA / "gg_adult.h5ad"
    a.write_h5ad(out, compression="gzip")
    print(f"wrote {out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--label_col", nargs="+", required=True,
                   help="One or more obs columns; multiple are joined with ' | '.")
    p.add_argument("--min_label_cells", type=int, default=10)
    p.add_argument("--drop_labels", nargs="*", default=[])
    main(**vars(p.parse_args()))
