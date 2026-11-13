#!/usr/bin/env python
"""Compare the visual phyloP assessment of MAFB-locus DARs with programmatic scores.

Input: results/manual_phylop/manual_phylop.tsv -- the "manual phylop" tab of the "DARs of
interest" sheet: 58 SnapATAC DARs at MAFB (34 RA, 24 G4) rated neutral / weak / moderate / strong
for whole-tree and oscine constraint by eye in the browser. NC_042585.1 = chr20 (lonStrDom2).
The G4 coordinates are approximate (the G4 DAR track was misbehaving).

Each region is scored over its full extent -- what the browser shows -- on every phyloP track:
mean phyloP, % of bases with phyloP > 1 and > 2, and the mean of its best 50-bp window (the
visual eye is drawn to a conserved block). Reported:
  * Spearman correlation of each score with the visual rating (0-3)
  * RA vs G4 on these same regions, by eye (share rated moderate/strong) and by score (MW)
  * overlap of each region with the loci_focus DAR elements
Writes results/manual_phylop/{manual_phylop_scored.tsv, agreement.tsv}.
"""
import re, sys
import numpy as np, pandas as pd, pyBigWig
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
O = f"{F}/results/manual_phylop"
T = "/private/warm-archive/colquittlab/assembly/tracks_lark/phyloP"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
TRACKS = {"whole tree (matched)": f"{PB}/wholetree_perbase_lsd2.bw",
          "oscine lineage Anc222 (matched)": f"{PB}/subtree222_perbase_lsd2.bw",
          "avian all (3-model)": f"{T}/phyloP_avian363_3model_all.bw",
          "vert. all (3-model)": f"{T}/phyloP_mamAvi605_3model_all.bw",
          "oscine Anc218 (3-model)": f"{T}/phyloP_avian363_3model_birdAnc218.bw"}
RATE = {"neutral": 0, "weak": 1, "moderate": 2, "modeate": 2, "strong": 3}
WHOLE = ["whole tree (matched)", "avian all (3-model)", "vert. all (3-model)"]
OSC = ["oscine lineage Anc222 (matched)", "oscine Anc218 (3-model)"]


def parse(r):
    c, a, b = re.match(r"\s*(\S+):([\d,]+)-([\d,]+)", r).groups()
    return ("chr20" if c == "NC_042585.1" else c), int(a.replace(",", "")) - 1, int(b.replace(",", ""))


def main():
    m = pd.read_csv(f"{O}/manual_phylop.tsv", sep="\t", dtype=str).fillna("")
    m[["chrom", "start", "end"]] = [parse(r) for r in m.region]
    m["whole_r"] = m.whole.str.strip().map(RATE); m["oscine_r"] = m.oscine.str.strip().map(RATE)
    bws = {k: pyBigWig.open(p) for k, p in TRACKS.items()}
    for k, bw in bws.items():
        cols = {"mean": [], "pct1": [], "pct2": [], "best50": []}
        for c, s, e in zip(m.chrom, m.start, m.end):
            v = np.nan_to_num(np.array(bw.values(c, s, e), float))
            cols["mean"].append(v.mean()); cols["pct1"].append(100 * (v > 1).mean()); cols["pct2"].append(100 * (v > 2).mean())
            w = np.convolve(v, np.ones(50) / 50, "valid") if len(v) >= 50 else v
            cols["best50"].append(w.max())
        for s_, x in cols.items(): m[f"{k}|{s_}"] = x
    el = pd.read_csv(f"{F}/allelem_cls.bed", sep="\t", header=None, names=["c", "s", "e", "n", "k"])
    el = el[el.c == "chr20"]
    m["loci_focus_DARs"] = [",".join(f"{n}({k[:2]})" for n, k in zip(*[el[(el.e > s) & (el.s < e)][x] for x in ("n", "k")]))
                            for s, e in zip(m.start, m.end)]
    m["width"] = m.end - m.start
    m.to_csv(f"{O}/manual_phylop_scored.tsv", sep="\t", index=False, float_format="%.3f")

    rows = []
    for rating, tracks in (("whole_r", list(TRACKS)), ("oscine_r", list(TRACKS))):
        for k in tracks:
            for s_ in ("mean", "pct1", "pct2", "best50"):
                rho, p = stats.spearmanr(m[rating], m[f"{k}|{s_}"])
                rows.append(dict(visual=rating.replace("_r", ""), track=k, score=s_, spearman=rho, p=p))
    A = pd.DataFrame(rows); A.to_csv(f"{O}/agreement.tsv", sep="\t", index=False, float_format="%.4g")
    print("== agreement of each score with the visual rating (Spearman rho, best scores per track)")
    for vis in ("whole", "oscine"):
        print(f"  visual '{vis}':")
        for k in TRACKS:
            a = A[(A.visual == vis) & (A.track == k)].set_index("score")
            print(f"    {k:32s} " + "  ".join(f"{s_} {a.loc[s_, 'spearman']:+.2f} (p={a.loc[s_, 'p']:.1g})" for s_ in ("mean", "pct1", "pct2", "best50")))

    print("\n== RA vs G4 on these 58 regions")
    for vis in ("whole_r", "oscine_r"):
        for k in ("RA", "G4"):
            x = m[m["class"] == k][vis]
            print(f"  visual {vis[:-2]:6s} {k}: n={len(x)}  rating mean {x.mean():.2f}  moderate/strong {100*(x>=2).mean():.0f}%  strong {100*(x==3).mean():.0f}%")
        a, b = m[m["class"] == "RA"][vis], m[m["class"] == "G4"][vis]
        print(f"    MW on ratings G4 vs RA: p = {stats.mannwhitneyu(b, a).pvalue:.3g}")
    for k in TRACKS:
        for s_ in ("pct1", "best50"):
            a, b = m[m["class"] == "RA"][f"{k}|{s_}"], m[m["class"] == "G4"][f"{k}|{s_}"]
            print(f"  {k:32s} {s_:6s} RA median {a.median():6.2f}  G4 median {b.median():6.2f}  MW p = {stats.mannwhitneyu(b, a).pvalue:.3g}")
    print(f"\n  regions overlapping a loci_focus DAR element: {(m.loci_focus_DARs != '').sum()} of {len(m)}")
    print(f"  wrote {O}/manual_phylop_scored.tsv, agreement.tsv")


if __name__ == "__main__":
    main()
