"""Per-cell cross-species alignment from SAMap's own stitched manifold.

Replaces the SATURN-embedding conservation score, which failed validation: three
different summary statistics computed on SATURN's 256-d space (iLISI-based
conservation, reciprocal-MNN k=20, k=50) were all UNCORRELATED with SAMap's cluster
alignment (Spearman -0.07, -0.14, -0.14; all p > 0.28), and all three ranked Astro/OPC
near the bottom where SAMap puts them near the top. The disagreement is a property of
the embedding, not of the statistic, so the fix is to measure on SAMap's manifold
rather than to re-metricise SATURN's.

Definition taken verbatim from samap.analysis._compute_csim, so this is SAMap's own
quantity rather than a new one:

    score(cell) = (sum of that cell's CROSS-SPECIES edge weights in
                   samap.adata.obsp['connectivities']) / uns['mapping_K']

_compute_csim averages exactly this over the cells of a cluster to get the cluster
alignment score (before a max(CSIM, CSIM.T) symmetrisation, which has no cell-level
analogue -- so cluster means here match the UN-symmetrised CSIM row sums, and the
published cluster scores are recovered from the same numbers).

Also reported per cell:
    top_ref_cluster  the other-species cluster holding the most of that cell's
                     cross-species edge weight
    specificity      that cluster's share of the cell's total cross-species weight,
                     i.e. how concentrated the mapping is (1 = all weight on one type)
"""
import argparse
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
CMAP = "magma_r"
BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly")


def read_obs_col(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):                      # categorical
        cats = g["categories"][:]
        codes = g["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        return cats[codes]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


def main(joint_h5ad: Path, out_dir: Path, query_species: str, ref_species: str,
         query_h5ad: Path, native_umap_key: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(joint_h5ad, "r") as f:
        obs_names = read_obs_col(f, "_index")
        species = read_obs_col(f, "species")
        clusters = read_obs_col(f, "cluster;cluster_mapping_scores")   # holds sid_label
        K = float(np.array(f["uns"]["mapping_K"]).ravel()[0])
        c = f["obsp"]["connectivities"]
        C = sp.csr_matrix((c["data"][:], c["indices"][:], c["indptr"][:]),
                          shape=tuple(c.attrs["shape"]))
        umap = f["obsm"]["X_umap"][:]
    n = len(obs_names)
    print(f"joint {n:,} cells | mapping_K={K:.0f} | connectivities nnz {C.nnz:,}")

    qm = species == query_species
    rm = species == ref_species
    print(f"{query_species} {qm.sum():,} | {ref_species} {rm.sum():,}")

    # --- keep only cross-species edges, exactly as _compute_csim does ---
    Cc = C.tocoo()
    keep = species[Cc.row] != species[Cc.col]
    Xc = sp.csr_matrix((Cc.data[keep], (Cc.row[keep], Cc.col[keep])), shape=C.shape)
    print(f"cross-species edges: {Xc.nnz:,} of {C.nnz:,} ({Xc.nnz/C.nnz:.1%})")

    alignment = np.asarray(Xc.sum(axis=1)).ravel() / K

    # --- per-cell top target cluster and how concentrated the mapping is ---
    ref_cl = clusters.copy()
    cats, codes = np.unique(ref_cl, return_inverse=True)
    Xq = Xc[qm]                                        # query rows only
    col_code = codes[np.arange(n)]
    M = sp.csr_matrix((Xq.data, col_code[Xq.indices], Xq.indptr),
                      shape=(Xq.shape[0], cats.size))
    M.sum_duplicates()
    Md = M.toarray()
    tot = Md.sum(1)
    top_i = Md.argmax(1)
    top_lab = cats[top_i]
    top_w = Md[np.arange(Md.shape[0]), top_i]
    with np.errstate(invalid="ignore", divide="ignore"):
        specificity = np.where(tot > 0, top_w / np.maximum(tot, 1e-12), 0.0)

    strip = lambda s: s[len(ref_species) + 1:] if s.startswith(ref_species + "_") else s
    df = pd.DataFrame({
        "cell": obs_names[qm],
        "cluster": [strip(x) if False else (x[len(query_species)+1:] if x.startswith(query_species+"_") else x)
                    for x in clusters[qm]],
        "alignment": alignment[qm].round(4),
        "specificity": specificity.round(4),
        "top_ref_cluster": [strip(x) for x in top_lab],
        "n_cross_edges": np.diff(Xq.indptr),
    }).set_index("cell")
    df.to_csv(out_dir / "per_cell_samap_alignment.csv")

    summ = (df.groupby("cluster")[["alignment", "specificity"]]
              .agg(["mean", "median"]).round(4))
    summ.columns = ["_".join(c) for c in summ.columns]
    summ["n_cells"] = df.groupby("cluster").size()
    summ["frac_zero"] = df.groupby("cluster").alignment.apply(lambda s: (s == 0).mean()).round(3)
    summ = summ.sort_values("alignment_mean", ascending=False)
    summ.to_csv(out_dir / "per_cell_samap_by_cluster.csv")

    pd.set_option("display.width", 220)
    print(f"\nper-cell alignment: mean {alignment[qm].mean():.4f} "
          f"median {np.median(alignment[qm]):.4f} frac-zero {(alignment[qm]==0).mean():.3f}")
    print("\n=== most aligned clusters ===");  print(summ.head(12).to_string())
    print("\n=== least aligned clusters ==="); print(summ.tail(8).to_string())

    # --- sanity check: cluster means must reproduce SAMap's published scores ---
    pub = out_dir / "samap_finch_calls.csv"
    if pub.exists():
        from scipy.stats import spearmanr
        p = pd.read_csv(pub).set_index("finch_cluster")
        both = summ.join(p[["samap_score"]], how="inner")
        r, pv = spearmanr(both.alignment_mean, both.samap_score)
        print(f"\n[check] cluster-mean per-cell alignment vs published samap_score: "
              f"spearman {r:+.3f} (p={pv:.2g}, n={len(both)})")
        print("        (should be strongly positive -- same underlying quantity)")

    # ---------------- plots ----------------
    def paint(coords, values, title, path, cbar_label, vmin=None, vmax=None, s=2.0):
        fig, ax = plt.subplots(figsize=(6.4, 5.6))
        order = np.argsort(values)
        h = ax.scatter(coords[order, 0], coords[order, 1], c=values[order], s=s,
                       cmap=CMAP, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        for sd in ax.spines.values(): sd.set_visible(False)
        cb = fig.colorbar(h, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label(cbar_label, fontsize=8); cb.ax.tick_params(labelsize=7)
        fig.savefig(path, bbox_inches="tight")
        fig.savefig(str(path).replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    vmax = float(np.quantile(df.alignment, 0.99))       # long right tail
    U = umap[qm]
    paint(U, df.alignment.values, f"{query_species} cells, SAMap manifold — alignment",
          out_dir / "umap_samap_alignment.pdf", "SAMap alignment (cross-species edge wt / K)",
          vmin=0, vmax=vmax)
    paint(U, df.specificity.values, f"{query_species} cells, SAMap manifold — mapping specificity",
          out_dir / "umap_samap_specificity.pdf", "top-cluster share of cross-species weight",
          vmin=0, vmax=1)

    if query_h5ad.exists():
        import anndata as ad
        qa = ad.read_h5ad(query_h5ad, backed="r")
        if native_umap_key in qa.obsm:
            common = qa.obs_names.intersection(df.index)
            print(f"\nnative UMAP '{native_umap_key}': matched {len(common):,}/{qa.n_obs:,}")
            V = pd.DataFrame(np.asarray(qa.obsm[native_umap_key]), index=qa.obs_names).loc[common].values
            sub = df.loc[common]
            paint(V, sub.alignment.values, f"{query_species} developmental UMAP — SAMap alignment",
                  out_dir / "umap_native_samap_alignment.pdf",
                  "SAMap alignment (cross-species edge wt / K)", vmin=0, vmax=vmax, s=3.0)
            paint(V, sub.specificity.values, f"{query_species} developmental UMAP — mapping specificity",
                  out_dir / "umap_native_samap_specificity.pdf",
                  "top-cluster share of cross-species weight", vmin=0, vmax=1, s=3.0)
        else:
            print(f"native UMAP key '{native_umap_key}' not found; have {list(qa.obsm)}")
    print(f"\nwrote -> {out_dir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--joint_h5ad", type=Path, default=BASE / "results" / "samap_joint.h5ad")
    p.add_argument("--out_dir", type=Path, default=BASE / "results")
    p.add_argument("--query_species", default="bf")
    p.add_argument("--ref_species", default="gg")
    p.add_argument("--query_h5ad", type=Path,
                   default=Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data/bf_adult_glut.h5ad"))
    p.add_argument("--native_umap_key", default="X_umap_seurat")
    main(**vars(p.parse_args()))
