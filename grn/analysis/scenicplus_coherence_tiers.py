#!/usr/bin/env python3
"""Concordance tiers for every SCENIC+ eRegulon in every config. DESCRIPTIVE, not a quality grade.

Do not use these to discard regulators. Discordance can be real biology: a TF may act through cofactors, on poised
or primed chromatin, or with different roles in different lineages, so its expression, its regions' accessibility and
its targets' expression need not line up across cell types. A discordant eRegulon is one to look at more closely
(which cell types, and why), not a worse one. Stable AND concordant is simply the easiest kind to interpret.

An eRegulon's regions come from motif enrichment (+ region-to-gene links) and its target genes from TF-to-gene
expression links; nothing forces the two to be active in the same cells, so a regulon can be "region-active" in one
cell type and "gene-active" in another (EMX2: regions peak in GABA-LGE, genes in Astro/RA/NSC). Two cell-type-level
checks (mean AUC per cluster, 27 clusters), for every eRegulon:

    rho_RG  Spearman between the gene-based and region-based AUC profiles
            -- do the regulon's regions and its target genes agree on where it is active?
    rho_TF  Spearman between the TF's own expression profile and the gene-based AUC profile
            -- is the TF expressed where its target genes are active?

A check passes at rho >= 0.45, about one-sided p < 0.01 for 27 clusters. The clusters are related lineages, not
independent, so read the p-values as a screen, not a test. Tiers (positive, +/+, eRegulons only):

    Tier 1  both checks pass      Tier 2  one passes      Tier 3  neither

Negative (-/-) eRegulons get the metrics but no tier (the expected sign of the TF-target relation differs).

Per TF, across the correct-input sweep (configs 1-11; 12/13 are the legacy-input controls and 14/15 topic controls,
kept in the eRegulon table but not in the TF tiers), stability x coherence:

    A  present in >= 9 of 11 configs and Tier 1 in >= half of the configs where present   stable, concordant
    B  present in >= 9 of 11 but Tier 1 in < half                                          stable, discordant
    C  present in < 9 and Tier 1 in >= half                                                variable, concordant
    D  everything else                                                                     variable, discordant

    python scenicplus_coherence_tiers.py
"""
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/")
RES = B + "results/"
CONFIGS = [f"config{n}" for n in range(1, 14)]
SWEEP = [f"config{n}" for n in range(1, 12)]
RHO_PASS = 0.45
PRESENT_MIN = 9  # of the 11 sweep configs


def s(x):
    return x.decode() if isinstance(x, bytes) else x


# TF expression per cluster: new run's adata, .raw = counts -> normalize to 1e4 per cell -> log1p -> cluster mean
a = ad.read_h5ad(B + "anndata_rna/adata.h5ad")
raw = csr_matrix(a.raw.X)
tot = np.asarray(raw.sum(1)).ravel()
norm = raw.multiply(1e4 / tot[:, None]).tocsr()
norm.data = np.log1p(norm.data)
clusters = sorted(a.obs.cluster.unique())
onehot = csr_matrix((np.ones(len(a)), (np.arange(len(a)), [clusters.index(c) for c in a.obs.cluster])), shape=(len(a), len(clusters)))
expr = pd.DataFrame((onehot.T @ norm).toarray() / np.asarray(onehot.sum(0)).ravel()[:, None],
                    index=clusters, columns=list(a.raw.var_names))


def spear(x, y):
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan
    r = spearmanr(x, y, alternative="greater")
    return float(r.statistic), float(r.pvalue)


rows = []
for cfg in CONFIGS:
    full = f"{RES}{cfg}/outs/scplusmdata.h5mu"
    f = h5py.File(full if Path(full).exists() else full.replace("scplusmdata.h5mu", "scplusmdata_slim.h5mu"), "r")  # slim: same layout
    obs = f["mod"]["scRNA_counts"]["obs"]
    cell2cl = dict(zip((s(x) for x in obs["_index"][:]), np.array([s(c) for c in obs["cluster"]["categories"][:]])[obs["cluster"]["codes"][:]]))
    for kind in ("direct", "extended"):
        mats = {}
        for basis in ("gene", "region"):
            g = f["mod"][f"{kind}_{basis}_based_AUC"]
            cells = [s(x) for x in g["obs"]["Cell"][:]]
            cl = np.array([cell2cl[c] for c in cells])
            X = np.asarray(g["X"][:], dtype=np.float32)
            names = [s(x) for x in g["var"]["_index"][:]]
            M = np.vstack([X[cl == c].mean(0) for c in clusters])  # clusters x regulons
            mats[basis] = (pd.DataFrame(M, index=clusters, columns=names))
        gene, region = mats["gene"], mats["region"]
        key = lambda n: n.rsplit("_(", 1)[0]
        rkey = {key(n): n for n in region.columns}
        for n in gene.columns:
            k = key(n)
            if k not in rkey:
                continue
            tf = k.split("_")[0]
            sign = re.search(r"_(\+/\+|-/-|\+/-|-/\+)$", k).group(1)
            gp, rp = gene[n].values, region[rkey[k]].values
            rho_rg, p_rg = spear(gp, rp)
            rho_tf, p_tf = spear(expr[tf].reindex(clusters).values, gp) if tf in expr.columns else (np.nan, np.nan)
            tier = None
            if sign == "+/+":
                passes = int(rho_rg >= RHO_PASS) + int(rho_tf >= RHO_PASS) if not (np.isnan(rho_rg) or np.isnan(rho_tf)) else \
                    int(rho_rg >= RHO_PASS) + int(rho_tf >= RHO_PASS)
                tier = {2: 1, 1: 2, 0: 3}[passes]
            rows.append({
                "config": cfg, "TF": tf, "type": kind, "sign": sign,
                "n_target_genes": int(re.search(r"\((\d+)g\)", n).group(1)), "n_target_regions": int(re.search(r"\((\d+)r\)", rkey[k]).group(1)),
                "rho_gene_vs_region": round(rho_rg, 3), "rho_TF_expr_vs_gene": round(rho_tf, 3), "concordance_tier": tier,
                "gene_AUC_peak": clusters[int(np.argmax(gp))], "region_AUC_peak": clusters[int(np.argmax(rp))],
                "TF_expr_peak": clusters[int(np.argmax(expr[tf].reindex(clusters).values))] if tf in expr.columns else None,
            })
    f.close()
    print(cfg, "done", flush=True)

E = pd.DataFrame(rows)
E.to_csv(HERE / "scenicplus_coherence_eregulons.csv", index=False)

# per TF, across the sweep (positive eRegulons)
P = E[(E.sign == "+/+") & E.config.isin(SWEEP)]
out = []
for tf, g in P.groupby("TF"):
    per_cfg = g.groupby("config").concordance_tier.min()                     # best tier of any +/+ eRegulon of this TF in a config
    present = per_cfg.size
    tier1_cfg = int((per_cfg == 1).sum())
    direct_cfg = int(g[g.type == "direct"].config.nunique())
    cls = ("A" if present >= PRESENT_MIN and tier1_cfg >= present / 2 else
           "B" if present >= PRESENT_MIN else
           "C" if tier1_cfg >= present / 2 else "D")
    out.append({"TF": tf, "stability_concordance_class": cls, "n_configs_present": present, "n_configs_tier1": tier1_cfg,
                "n_configs_with_direct": direct_cfg, "median_rho_gene_vs_region": round(g.rho_gene_vs_region.median(), 3),
                "median_rho_TF_expr_vs_gene": round(g.rho_TF_expr_vs_gene.median(), 3),
                "TF_expr_peak": g.TF_expr_peak.mode().iat[0] if g.TF_expr_peak.notna().any() else None,
                "gene_AUC_peak": g.gene_AUC_peak.mode().iat[0], "region_AUC_peak": g.region_AUC_peak.mode().iat[0]})
T = pd.DataFrame(out).sort_values(["stability_concordance_class", "n_configs_tier1"], ascending=[True, False])
T.to_csv(HERE / "scenicplus_coherence_tf_tiers.csv", index=False)

print("\neRegulon tiers (+/+, sweep configs 1-11):", P.concordance_tier.value_counts().sort_index().to_dict())
print("TF classes:", T["stability_concordance_class"].value_counts().sort_index().to_dict())
for k in "ABCD":
    sub = T[T["stability_concordance_class"] == k]
    print(f"  {k}: {len(sub)} TFs", ("; " + ", ".join(sub.TF.head(30))) if k == "A" else "")
for tf in ("EMX2", "LHX2", "MAFB"):
    r = T[T.TF == tf]
    print(tf, r[["stability_concordance_class", "n_configs_present", "n_configs_tier1", "median_rho_gene_vs_region", "median_rho_TF_expr_vs_gene"]].to_dict("records") or "no +/+ eRegulon in the sweep")
