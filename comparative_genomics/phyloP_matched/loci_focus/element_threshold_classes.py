#!/usr/bin/env python
"""Threshold classes on whole-tree and oscine mean phyloP, compared with the Ward groups.

Features are those of plot_element_cluster.R noacc: each lobe's whole-tree and oscine-lineage mean
phyloP over its full extent, accessibility-weighted (whole_ext_mean, osc_ext_mean). Thresholds are
set on the visual ratings as in element_rule_classes.py (rated moderate/strong vs neutral/weak;
Youden's J; best-scoring overlapping lobe per rated region) and are absolute, so a depth filter only
removes lobes and never moves the cut. Classes: deep+oscine, deep (split at an avian mean of 3.0 into
strong avian and moderate avian -- the top Ward split of the deep lobes, which falls in an empty gap
between 2.63 and 3.21), oscine, weak. Each is cross-
tabulated against the Ward k=3 groups (element_cluster_k3_noacc[_depth<T>].tsv) without and with the
depth filters (every 50-bp window of the lobe aligned in >= T of oscine and non-oscine genomes).

    python element_threshold_classes.py
Writes results/mafb_domain/element_threshold_classes.tsv.
"""
import os
import numpy as np, pandas as pd
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
R = f"{F}/results/mafb_domain"
SC = {"whole": "whole_ext_mean", "oscine": "osc_ext_mean"}
ORDER = ["deep+oscine", "strong avian", "moderate avian", "oscine", "weak"]
STRONG = 3.0                       # avian (whole-tree) mean splitting deep-only lobes: the Ward split, in a gap 2.63-3.21
UP = {"moderate", "strong"}


def youden(score, up):
    best = (-1, None)
    for t in np.unique(np.round(score, 2)):
        call = score >= t; j = call[up].mean() + (~call[~up]).mean() - 1
        if j > best[0]: best = (j, t)
    return best


def main():
    S = pd.read_csv(f"{R}/element_rule_classes.tsv", sep="\t")
    m = pd.read_csv(f"{F}/results/manual_phylop/manual_phylop.tsv", sep="\t", dtype=str).fillna("")
    m = m[m.whole.isin(["neutral", "weak", "moderate", "strong"])]; rows = []
    for r in m.itertuples(index=False):
        a, b = r.region.strip().split(":")[-1].replace(",", "").split("-"); a, b = int(a) - 1, int(b)
        h = S[(S.ext_start < b) & (S.ext_end > a) & S.cls.isin([r[3], "both"])]
        if len(h): rows.append(dict(whole_r=r.whole, oscine_r=r.oscine, **{t: h[c].max() for t, c in SC.items()}))
    M = pd.DataFrame(rows); T = {}
    for t in SC:
        up = M[f"{t}_r"].isin(UP).values; j, thr = youden(M[t].values, up); T[t] = thr; call = M[t].values >= thr
        print(f"  {t:6s} mean phyloP threshold {thr:+.2f}: agreement {(call == up).mean():.0%} (sens {call[up].mean():.0%}, spec {(~call[~up]).mean():.0%}), J {j:.2f}")
    d, o = S[SC["whole"]] >= T["whole"], S[SC["oscine"]] >= T["oscine"]
    s3 = S[SC["whole"]] >= STRONG
    S["thr_class"] = np.select([d & o, d & s3, d, o], ORDER[:4], ORDER[4])
    S[["element", "cls", "chrom", "ext_start", "ext_end", "whole_ext_mean", "osc_ext_mean", "osc_depth_min50", "nonosc_depth_min50", "thr_class"]].assign(
        thr_whole=T["whole"], thr_oscine=T["oscine"]).to_csv(f"{R}/element_threshold_classes.tsv", sep="\t", index=False, float_format="%.3f")
    for tag, dmin in (("", None), ("_depth0.8", 0.8), ("_depth0.9", 0.9)):
        x = S if dmin is None else S[(S.osc_depth_min50 >= dmin) & (S.nonosc_depth_min50 >= dmin)]
        f = f"{R}/element_cluster_k3_noacc{tag}.tsv"
        print(f"\n== {'no depth filter' if dmin is None else f'depth >= {dmin}'}: {len(x)} lobes")
        print(pd.crosstab(pd.Categorical(x.thr_class, ORDER), x.cls, rownames=["threshold"], margins=True).to_string())
        if os.path.exists(f):
            w = x.merge(pd.read_csv(f, sep="\t")[["element", "group"]], on="element")
            print("  threshold (rows) vs Ward k=3 (cols):")
            print(pd.crosstab(pd.Categorical(w.thr_class, ORDER), w.group, rownames=["threshold"]).to_string())
            tc = w.thr_class.replace({"strong avian": "deep", "moderate avian": "deep"})
            agree = ((tc == w.group) | ((tc == "deep+oscine") & w.group.isin(["deep", "oscine"]))).mean()
            print(f"  same call (deep+oscine counted as matching deep or oscine): {agree:.0%}")
        t2 = x[x.cls.isin(["RA", "G4"])]
        print(f"  RA vs G4 across threshold classes: chi-square p = {stats.chi2_contingency(pd.crosstab(t2.thr_class, t2.cls))[1]:.2g}")


if __name__ == "__main__":
    main()
