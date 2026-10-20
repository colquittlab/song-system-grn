#!/usr/bin/env python3
"""Summary of tf2g_multiseed_stability.py (AR vs NR3C2, RA vs C1H-1 all runs, config2): which TF-to-gene links are stable across random seeds?

A link "AR -> gene" is called STABLE if AR is among the top-10 predictors of that gene in at least STABLE_MIN of the N seeds. Reported:
  1. seed-to-seed consistency of AR's and NR3C2's importance across genes
  2. numbers of stable links, overlap of the two TFs' stable links, and how many of each TF's eRegulon targets are stable
  3. MAFB: AR's and NR3C2's rank / importance across seeds, and the TFs that are stable top predictors of MAFB
  4. drop-one: with the seed noise averaged out, does AR gain importance when NR3C2 is removed (and the reverse)?
Writes tf2g_stability_mafb_predictors.csv and tf2g_stability_summary.csv next to this script.

    python tf2g_stability_summary.py
"""
import itertools
import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
        "glut-ra_glut-c1h-1_allruns_hybrid/results/config2/")
S = ROOT + "tf2g_stability/"
TOPK, STABLE_MIN = 10, 8
R = pd.read_csv(S + "tf2g_ar_nr3c2_by_gene_seed.csv")
M = pd.read_csv(S + "tf2g_mafb_predictors_by_seed.csv")
rho = pd.read_csv(S + "rho_ar_nr3c2_by_gene.csv", index_col=0)
er = pd.concat([pd.read_table(ROOT + "outs/" + f, usecols=["TF", "Gene"]) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
n_seeds = R.seed.nunique()
full = R[R.tag == "full"].copy()
for t in ("AR", "NR3C2"):
    full[f"{t}_share"] = full[f"{t}_imp"] / full.total
    full[f"{t}_top"] = (full[f"{t}_rank"] <= TOPK).astype(int)
out = []


def rec(section, name, value):
    out.append(dict(section=section, metric=name, value=value)); return value


print(f"{n_seeds} seeds, {full.gene.nunique()} genes, top-{TOPK} in >= {STABLE_MIN} seeds = stable\n")
print("1. Seed-to-seed consistency of the TF's importance across genes (mean pairwise Spearman over seed pairs):")
for t in ("AR", "NR3C2"):
    P = full.pivot(index="gene", columns="seed", values=f"{t}_imp").fillna(0)
    rs = [spearmanr(P[a], P[b])[0] for a, b in itertools.combinations(P.columns, 2)]
    sd = P.std(axis=1); mu = P.mean(axis=1)
    print(f"   {t:6} mean pairwise rho {np.mean(rs):.3f}; genes where {t} is in the top-{TOPK} in >=1 seed: {int((P.index.isin(full[full[f'{t}_top']==1].gene)).sum())}")
    rec("consistency", f"{t}_mean_pairwise_spearman", round(float(np.mean(rs)), 3))

freq = {t: full.groupby("gene")[f"{t}_top"].sum() for t in ("AR", "NR3C2")}
mean_share = {t: full.groupby("gene")[f"{t}_share"].mean() for t in ("AR", "NR3C2")}
stable = {t: set(freq[t][freq[t] >= STABLE_MIN].index) for t in freq}
stable_pos = {t: {g for g in stable[t] if rho.loc[g, t] > 0.05} for t in stable}
print("\n2. Stable links (all genes):")
for t in ("AR", "NR3C2"):
    print(f"   {t:6} {len(stable[t]):4d} stable links ({len(stable_pos[t])} with Spearman rho > 0.05); top-{TOPK} in >=1 seed for {int((freq[t] >= 1).sum())} genes")
    rec("stable_links", f"{t}_n_stable", len(stable[t])); rec("stable_links", f"{t}_n_stable_rho_pos", len(stable_pos[t]))
shared = stable["AR"] & stable["NR3C2"]
print(f"   genes stable for BOTH: {len(shared)} (of {len(stable['AR'])} AR, {len(stable['NR3C2'])} NR3C2)")
rec("stable_links", "n_stable_both", len(shared))
for t in ("AR", "NR3C2"):
    tg = {g for g in set(er[er.TF == t].Gene) if g in set(full.gene) and g != t}
    st = len(tg & stable[t]); s1 = int((freq[t].reindex(list(tg)).fillna(0) >= 1).sum())
    print(f"   {t:6} eRegulon targets {len(tg)}: stable (top-{TOPK} in >={STABLE_MIN}/{n_seeds} seeds) {st} ({100*st/len(tg):.0f}%); top-{TOPK} in >=1 seed {s1} ({100*s1/len(tg):.0f}%)")
    rec("eregulon_targets", f"{t}_targets", len(tg)); rec("eregulon_targets", f"{t}_targets_stable", st); rec("eregulon_targets", f"{t}_targets_top10_in_ge1_seed", s1)

print("\n3. MAFB as a target:")
mf = M[M.tag == "full"]
tab = (mf.groupby("TF").agg(mean_importance=("importance", "mean"), sd=("importance", "std"), mean_rank=("rank", "mean"),
       top10_seeds=("rank", lambda r: int((r <= TOPK).sum())), seeds_present=("seed", "nunique")).sort_values("mean_importance", ascending=False))
tab["top10_seeds"] = tab.top10_seeds.astype(int)
for t in ("AR", "NR3C2"):
    r = tab.loc[t] if t in tab.index else None
    if r is None:
        print(f"   {t}: never selected"); continue
    print(f"   {t:6} importance {r.mean_importance:.3f} +/- {r.sd:.3f}, mean rank {r.mean_rank:.1f} of {int(mf.groupby('seed').TF.nunique().mean())}, top-{TOPK} in {r.top10_seeds}/{n_seeds} seeds, present in {int(r.seeds_present)} seeds")
    rec("MAFB", f"{t}_top10_seeds", int(r.top10_seeds)); rec("MAFB", f"{t}_mean_rank", round(float(r.mean_rank), 1)); rec("MAFB", f"{t}_mean_importance", round(float(r.mean_importance), 3))
inputs = set(er[er.Gene == "MAFB"].TF) - {"MAFB"}
tab["rho_with_MAFB_in_eRegulon_table"] = [np.nan] * len(tab)
tab["in_config2_eRegulon_inputs"] = tab.index.isin(inputs)
tab.round(3).head(25).to_csv(f"{HERE}/tf2g_stability_mafb_predictors.csv")
print(f"   TFs that are stable top-{TOPK} predictors of MAFB (top-{TOPK} in >= {STABLE_MIN} seeds): {list(tab[tab.top10_seeds >= STABLE_MIN].index)}")
print(f"   MAFB regulators in the config2 eRegulon tables: {sorted(inputs)}")
print("   top 12 by mean importance:")
print(tab.head(12)[["mean_importance", "sd", "mean_rank", "top10_seeds", "in_config2_eRegulon_inputs"]].round(3).to_string())

print("\n4. Drop-one, seed noise averaged out (genes = AR targets + NR3C2 targets + MAFB):")
sub = R[R.gene.isin(set(R[R.tag == "noNR3C2"].gene))].copy()
for t in ("AR", "NR3C2"):
    sub[f"{t}_share"] = sub[f"{t}_imp"] / sub.total
for t, drop, tag in (("AR", "NR3C2", "noNR3C2"), ("NR3C2", "AR", "noAR")):
    a = sub[sub.tag == "full"].groupby("gene")[f"{t}_share"].agg(["mean", "std"])
    b = sub[sub.tag == tag].groupby("gene")[f"{t}_share"].agg(["mean", "std"])
    j = a.join(b, lsuffix="_full", rsuffix="_drop").dropna()
    se = np.sqrt(j.std_full ** 2 / n_seeds + j.std_drop ** 2 / n_seeds).replace(0, np.nan)
    z = (j.mean_drop - j.mean_full) / se
    print(f"   {t:6} without {drop:6}: mean share {j.mean_full.mean():.4f} -> {j.mean_drop.mean():.4f} ({100*(j.mean_drop.mean()/j.mean_full.mean()-1):+.1f}%); "
          f"genes with a significant gain (z>2) {int((z > 2).sum())}, loss (z<-2) {int((z < -2).sum())} of {len(j)}; median per-gene change {100*float(((j.mean_drop/j.mean_full.replace(0,np.nan))-1).median()):+.1f}%")
    rec("drop_one", f"{t}_without_{drop}_mean_share_change_pct", round(100 * float(j.mean_drop.mean() / j.mean_full.mean() - 1), 1))
    rec("drop_one", f"{t}_without_{drop}_genes_gain_z2", int((z > 2).sum())); rec("drop_one", f"{t}_without_{drop}_genes_loss_z2", int((z < -2).sum()))
    if "MAFB" in j.index and t == "AR":
        print(f"      MAFB: AR share {j.loc['MAFB','mean_full']:.4f} (sd {j.loc['MAFB','std_full']:.4f}) -> {j.loc['MAFB','mean_drop']:.4f} without NR3C2")
pd.DataFrame(out).to_csv(f"{HERE}/tf2g_stability_summary.csv", index=False)
