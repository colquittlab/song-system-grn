#!/usr/bin/env python3
"""Summary of the AR / NR3C2 in-silico perturbation (glut-ra_glut-c1h-1_allruns_hybrid, config2; perturbation_AR-NR3C2/).

Reads the cached simulations (perturb_<TF>_sv0|sv10_<iter>.parquet), the regressors and the eRegulon tables, and reports for each TF and mode
(KD = set to 0, OE = set to 10 counts; iteration 5 vs the original matrix, raw counts):
  - the predicted log2FC of MAFB and of the TF's own target genes, per cluster
  - how similar the AR and NR3C2 effects are (Spearman over genes; overlap of the most-affected genes)
  - the direction of the PCA flow field in each cluster relative to the RA -> C1H-1 axis (from the cached arrow fields, seed 0)
Writes ar_nr3c2_perturbation_summary.csv next to this script.

    python ar_nr3c2_perturbation_summary.py
"""
import os
import pickle
import warnings

import anndata
import mudata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
        "glut-ra_glut-c1h-1_allruns_hybrid/results/config2/")
P = ROOT + "perturbation_AR-NR3C2/"
TFS, MODES = ("AR", "NR3C2"), (("KD", "sv0"), ("OE", "sv10"))
RA, C1 = "Glut-CACNA1H-RA", "Glut-CACNA1H-1"

er = pd.concat([pd.read_table(ROOT + "outs/" + f, usecols=["TF", "Gene"]) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
reg = pickle.load(open(P + "regressors.pkl", "rb"))
h5 = ROOT + "outs/scplusmdata_perturb.h5mu"
rna = mudata.read(h5 + "/scRNA_counts")
cl = rna.obs["cluster"].values
targets = {tf: sorted(g for g in set(er[er.TF == tf].Gene) if g in reg and g != tf and "LOC" not in g) for tf in TFS}


def sim(tf, tag, i):
    return pd.read_parquet(f"{P}perturb_{tf}_{tag}_{i}.parquet", engine="fastparquet")


def l2fc(tf, tag):
    """per-gene, per-cluster log2FC of the cluster means (iteration 5 / original)"""
    p0, p5 = sim(tf, tag, 0), sim(tf, tag, 5)
    out = {}
    for c in (RA, C1):
        idx = np.where(cl == c)[0]
        b, a = p0.iloc[idx].mean(), p5.iloc[idx].mean()
        ok = b > 0.02
        out[c] = np.log2(a[ok] / b[ok])
    return pd.DataFrame(out)


rows, eff = [], {}
for tf in TFS:
    for mode, tag in MODES:
        L = l2fc(tf, tag)
        eff[(tf, mode)] = L
        own = L.reindex(targets[tf])
        for c in (RA, C1):
            rows.append(dict(TF=tf, mode=mode, cluster=c, n_targets=len(targets[tf]), n_targets_scored=int(own[c].notna().sum()),
                             MAFB_log2FC=round(float(L[c].get("MAFB", np.nan)), 3), own_targets_median_log2FC=round(float(own[c].median()), 4),
                             own_targets_pct_changed_gt0p1=round(100 * float((own[c].abs() > 0.1).mean()), 1),
                             own_targets_max_abs=round(float(own[c].abs().max()), 3)))
S = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print(S.to_string(index=False))

print("\nTop affected target genes of each TF (RA cluster, signed log2FC):")
for tf in TFS:
    for mode, _ in MODES:
        s = eff[(tf, mode)].reindex(targets[tf])[RA].dropna()
        top = s.reindex(s.abs().sort_values(ascending=False).index)[:6]
        print(f"  {tf} {mode}: " + ", ".join(f"{g} {v:+.2f}" for g, v in top.items()))

print("\nSimilarity of the AR and NR3C2 effects (RA cluster):")
for mode, _ in MODES:
    a, b = eff[("AR", mode)][RA], eff[("NR3C2", mode)][RA]
    both = a.index.intersection(b.index)
    both = both[np.isfinite(a[both].to_numpy()) & np.isfinite(b[both].to_numpy())]   # genes with a zero baseline give inf/nan
    rho = spearmanr(a[both], b[both])[0]
    ta, tb = set(a.abs().sort_values(ascending=False).index[:100]), set(b.abs().sort_values(ascending=False).index[:100])
    shared_t = sorted(set(targets["AR"]) & set(targets["NR3C2"]))
    print(f"  {mode}: Spearman over {len(both)} genes {rho:+.2f}; top-100 most-affected overlap {len(ta & tb)}; "
          f"AR and NR3C2 share {len(shared_t)} targets (of {len(targets['AR'])} / {len(targets['NR3C2'])})")
    S.attrs[f"rho_{mode}"] = rho

# flow-field direction relative to the RA -> C1H-1 axis (cached fields: PCA+Harmony of the AUC matrix, seed 0)
d_ = mudata.read(h5 + "/direct_gene_based_AUC"); e_ = mudata.read(h5 + "/extended_gene_based_AUC")
auc = anndata.concat([d_, e_], axis=1); auc.obs = rna.obs
sc.pp.pca(auc); sc.external.pp.harmony_integrate(auc, key="replicate")
emb = auc.obsm["X_pca_harmony"][:, :2]
ra_c, c1_c = emb[cl == RA].mean(0), emb[cl == C1].mean(0)
axis = (c1_c - ra_c) / np.linalg.norm(c1_c - ra_c)
width = np.ptp(emb[:, 0])
print(f"\nFlow-field direction (cosine of the mean arrow with the RA -> C1H-1 axis; arrows >= 2% of the longest, grid points within 20% of the centroid):")
arrow_rows = []
for tf in TFS:
    for mode, _ in (("KD", 0), ("OE", 0)):
        pk = f"{P}arrow_field_{mode.lower()}_{tf}_seed0.pkl"
        if not os.path.exists(pk):
            continue
        f = pickle.load(open(pk, "rb"))["field"]
        keep = f["mask"] & (f["scaled"] >= 0.02) & (f["dist"] > 0)
        for name, cen in (("RA", ra_c), ("C1H-1", c1_c)):
            near = keep & (np.linalg.norm(f["grid_xy"] - cen, axis=1) < 0.20 * width)
            if near.sum() == 0:
                print(f"  {tf} {mode} near {name}: no arrows"); continue
            m = f["uv"][near].mean(0)
            cos = float(np.dot(m, axis) / (np.linalg.norm(m) + 1e-12))
            arrow_rows.append(dict(TF=tf, mode=mode, near=name, n_arrows=int(near.sum()), cos_with_RA_to_C1H1_axis=round(cos, 2)))
            print(f"  {tf} {mode} near {name}: {int(near.sum())} arrows, cosine {cos:+.2f}")
out = pd.concat([S.assign(kind="effect"), pd.DataFrame(arrow_rows).assign(kind="flow")], ignore_index=True)
out.to_csv(f"{HERE}/ar_nr3c2_perturbation_summary.csv", index=False)
