#!/usr/bin/env python
"""Rule-based constraint classes for the MAFB-domain accessible elements, calibrated on the visual ratings.

Each lobe (element_profiles.py) is scored over its full extent, weighted by accessibility:
  whole   whole_ext_block   % of bases in 11-bp blocks with mean whole-tree phyloP > 1
  oscine  osc_ext_pct1      % of bases with oscine-lineage (Anc222) phyloP > 1
These were the scores that tracked the visual ratings best (manual_phylop_compare.py). Thresholds
are set where they best separate rated regions called moderate or strong from those called
neutral or weak (Youden's J; each rated region takes its best-scoring overlapping lobe, same
class or "both"). Classes:
  deep+oscine   whole >= Tw and oscine >= To
  deep          whole >= Tw
  oscine        oscine >= To
  weak          neither

    python element_rule_classes.py
Writes results/mafb_domain/element_rule_classes.tsv and prints thresholds, agreement and counts.
"""
import numpy as np, pandas as pd
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
R = f"{F}/results/mafb_domain"
SC = {"whole": "whole_ext_block", "oscine": "osc_ext_pct1"}
UP = {"moderate", "strong"}
ORDER = ["deep+oscine", "deep", "oscine", "weak"]


def coords(reg):
    a, b = reg.strip().split(":")[-1].replace(",", "").split("-")
    return int(a) - 1, int(b)


def youden(score, up):
    best = (-1, None)
    for t in np.unique(np.round(score, 1)):
        call = score >= t
        j = call[up].mean() + (~call[~up]).mean() - 1
        if j > best[0]: best = (j, t)
    return best


def main():
    S = pd.read_csv(f"{R}/element_summaries.tsv", sep="\t")
    m = pd.read_csv(f"{F}/results/manual_phylop/manual_phylop.tsv", sep="\t", dtype=str).fillna("")
    m = m[m.whole.isin(["neutral", "weak", "moderate", "strong"])]
    rows = []
    for r in m.itertuples(index=False):
        a, b = coords(r.region)
        h = S[(S.ext_start < b) & (S.ext_end > a) & S.cls.isin([r[3], "both"])]
        if h.empty: continue
        rows.append(dict(region=r.region.strip(), cls=r[3], whole_r=r.whole, oscine_r=r.oscine,
                         **{t: h[c].max() for t, c in SC.items()}))
    M = pd.DataFrame(rows)
    print(f"  {len(M)} of {len(m)} rated regions matched to a scored lobe")
    T = {}
    for t, c in SC.items():
        up = M[f"{t}_r"].isin(UP).values; sc = M[t].values
        j, thr = youden(sc, up); T[t] = thr; call = sc >= thr
        rank = M[f"{t}_r"].map({"neutral": 0, "weak": 1, "moderate": 2, "strong": 3})
        print(f"  {t:6s} {c:16s} threshold {thr:5.1f}   J {j:.2f}   agreement {(call == up).mean():.0%} "
              f"(sens {call[up].mean():.0%}, spec {(~call[~up]).mean():.0%})   Spearman with rating {stats.spearmanr(sc, rank)[0]:+.2f}")
        for x in M[call != up].itertuples():
            print(f"      disagree: {x.region:34s} {x.cls:2s} rated {getattr(x, t + '_r'):8s} score {getattr(x, t):5.1f}")
    S["deep"] = S[SC["whole"]] >= T["whole"]; S["osc"] = S[SC["oscine"]] >= T["oscine"]
    S["rule"] = np.select([S.deep & S.osc, S.deep, S.osc], ORDER[:3], ORDER[3])
    S["thr_whole"], S["thr_oscine"] = T["whole"], T["oscine"]
    S.drop(columns=["deep", "osc"]).to_csv(f"{R}/element_rule_classes.tsv", sep="\t", index=False, float_format="%.3f")
    tab = pd.crosstab(S.rule, S.cls).reindex(ORDER).fillna(0).astype(int)
    print(tab.assign(total=tab.sum(1)).to_string())
    print(f"  class x rule Fisher-Freeman-Halton (RA vs G4, chi2 approx) p = {stats.chi2_contingency(tab[['RA', 'G4']])[1]:.3g}")


if __name__ == "__main__":
    main()
