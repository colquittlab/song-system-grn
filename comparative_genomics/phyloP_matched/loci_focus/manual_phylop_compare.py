#!/usr/bin/env python
"""Compare the visual phyloP assessment of MAFB-locus DARs with programmatic scores.

Input: results/manual_phylop/manual_phylop.tsv -- the "manual phylop" tab of the "DARs of
interest" sheet: 58 SnapATAC DARs at MAFB (34 RA, 24 G4) rated neutral / weak / moderate / strong
for whole-tree and oscine constraint by eye in the browser. NC_042585.1 = chr20 (lonStrDom2).
The G4 coordinates are approximate (the G4 DAR track was misbehaving).

Each region is scored over its full extent -- what the browser shows -- on every phyloP track:
mean phyloP, % of bases with phyloP > 1 and > 2, and the mean of its best 50-bp window (the
visual eye is drawn to a conserved block). The eye judges conservation where the region is
accessible, and DAR calls can extend into neutral flanks or sit beside a conserved block, so
each region is also scored guided by its up cluster's fragment signal (Glut-RA for RA, GABA-4-1
for G4) within the region:
  summit     the 200 bp about the fragment-signal summit (21-bp smoothed)
  hidiff     bases with up-cluster minus reference-cluster fragment difference > 5 (the
             programmatic filter; Glut-RA - Glut-Arco-1, GABA-4-1 - GABA-1-1)
  accw       mean phyloP and share > 1 weighted by the up cluster's fragment signal
Continuity (the eye discounts scattered constrained bases):
  block      % of bases whose 11-bp rolling-mean phyloP > 1 (and accessibility-weighted, accw_block)
  clust      % of bases with phyloP > 1 lying in runs >= 10 bp (1-bp gaps bridged)
  runmax     longest such run (bp)
Reported:
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
FR = "/private/warm-archive/colquittlab/assembly/tracks_lark/fragment"
UP = {"RA": ("Glut-RA", "Glut-Arco-1"), "G4": ("GABA-4-1", "GABA-1-1")}
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
    fb = {k: pyBigWig.open(f"{FR}/{k}-fragment-bin1-smoothNone.bw") for k in ("Glut-RA", "Glut-Arco-1", "GABA-4-1", "GABA-1-1")}
    fr = []
    for c, s, e, k in zip(m.chrom, m.start, m.end, m["class"]):
        up, ref = (np.nan_to_num(np.array(fb[x].values(c, s, e), float)) for x in UP[k])
        sm = np.convolve(up, np.ones(21) / 21, "same"); j = int(np.argmax(sm))
        fr.append(dict(up=up, diff=np.convolve(up - ref, np.ones(11) / 11, "same"), summit=j))
    m["summit"] = [s + f["summit"] for s, f in zip(m.start, fr)]
    m["n_hidiff"] = [int((f["diff"] > 5).sum()) for f in fr]
    for k, bw in bws.items():
        cols = {"mean": [], "pct1": [], "pct2": [], "best50": []}
        for c, s, e in zip(m.chrom, m.start, m.end):
            v = np.nan_to_num(np.array(bw.values(c, s, e), float))
            cols["mean"].append(v.mean()); cols["pct1"].append(100 * (v > 1).mean()); cols["pct2"].append(100 * (v > 2).mean())
            w = np.convolve(v, np.ones(50) / 50, "valid") if len(v) >= 50 else v
            cols["best50"].append(w.max())
        for s_, x in cols.items(): m[f"{k}|{s_}"] = x
        g = {"summit_pct1": [], "summit_mean": [], "hidiff_pct1": [], "hidiff_mean": [], "accw_pct1": [], "accw_mean": []}
        for (c, s, e), f in zip(zip(m.chrom, m.start, m.end), fr):
            v = np.nan_to_num(np.array(bw.values(c, s, e), float))
            a, b = max(0, f["summit"] - 100), min(len(v), f["summit"] + 100); vs = v[a:b]
            g["summit_pct1"].append(100 * (vs > 1).mean()); g["summit_mean"].append(vs.mean())
            hm = f["diff"] > 5
            g["hidiff_pct1"].append(100 * (v[hm] > 1).mean() if hm.any() else np.nan)
            g["hidiff_mean"].append(v[hm].mean() if hm.any() else np.nan)
            w = f["up"]; ws = w.sum()
            g["accw_pct1"].append(100 * ((v > 1) * w).sum() / ws if ws else np.nan)
            g["accw_mean"].append((v * w).sum() / ws if ws else np.nan)
        for s_, x in g.items(): m[f"{k}|{s_}"] = x
        cb = {"block": [], "accw_block": [], "clust": [], "runmax": []}
        for (c, s, e), f in zip(zip(m.chrom, m.start, m.end), fr):
            v = np.nan_to_num(np.array(bw.values(c, s, e), float))
            sm = np.convolve(v, np.ones(11) / 11, "same"); blk = sm > 1
            cb["block"].append(100 * blk.mean())
            w = f["up"]; cb["accw_block"].append(100 * (blk * w).sum() / w.sum() if w.sum() else np.nan)
            hi = v > 1
            hi = hi | (np.r_[False, hi[:-1]] & np.r_[hi[1:], False])        # bridge single-base gaps
            runs, n0 = [], 0
            for x in np.r_[hi, False]:
                if x: n0 += 1
                elif n0: runs.append(n0); n0 = 0
            cb["clust"].append(100 * sum(r for r in runs if r >= 10) / len(v))
            cb["runmax"].append(max(runs) if runs else 0)
        for s_, x in cb.items(): m[f"{k}|{s_}"] = x
    el = pd.read_csv(f"{F}/allelem_cls.bed", sep="\t", header=None, names=["c", "s", "e", "n", "k"])
    el = el[el.c == "chr20"]
    m["loci_focus_DARs"] = [",".join(f"{n}({k[:2]})" for n, k in zip(*[el[(el.e > s) & (el.s < e)][x] for x in ("n", "k")]))
                            for s, e in zip(m.start, m.end)]
    m["width"] = m.end - m.start
    m.to_csv(f"{O}/manual_phylop_scored.tsv", sep="\t", index=False, float_format="%.3f")

    rows = []
    for rating, tracks in (("whole_r", list(TRACKS)), ("oscine_r", list(TRACKS))):
        for k in tracks:
            for s_ in ("mean", "pct1", "pct2", "best50", "summit_pct1", "summit_mean", "hidiff_pct1", "hidiff_mean", "accw_pct1", "accw_mean",
                       "block", "accw_block", "clust", "runmax"):
                ok = m[f"{k}|{s_}"].notna()
                rho, p = stats.spearmanr(m[rating][ok], m[f"{k}|{s_}"][ok])
                rows.append(dict(visual=rating.replace("_r", ""), track=k, score=s_, spearman=rho, p=p))
    A = pd.DataFrame(rows); A.to_csv(f"{O}/agreement.tsv", sep="\t", index=False, float_format="%.4g")
    print("== agreement of each score with the visual rating (Spearman rho, best scores per track)")
    for vis in ("whole", "oscine"):
        print(f"  visual '{vis}':")
        for k in TRACKS:
            a = A[(A.visual == vis) & (A.track == k)].set_index("score")
            print(f"    {k:32s} " + "  ".join(f"{s_} {a.loc[s_, 'spearman']:+.2f}" for s_ in ("pct1", "summit_pct1", "accw_pct1", "block", "accw_block", "clust", "runmax")))

    print("\n== RA vs G4 on these 58 regions")
    for vis in ("whole_r", "oscine_r"):
        for k in ("RA", "G4"):
            x = m[m["class"] == k][vis]
            print(f"  visual {vis[:-2]:6s} {k}: n={len(x)}  rating mean {x.mean():.2f}  moderate/strong {100*(x>=2).mean():.0f}%  strong {100*(x==3).mean():.0f}%")
        a, b = m[m["class"] == "RA"][vis], m[m["class"] == "G4"][vis]
        print(f"    MW on ratings G4 vs RA: p = {stats.mannwhitneyu(b, a).pvalue:.3g}")
    for k in TRACKS:
        for s_ in ("pct1", "summit_pct1", "accw_pct1", "block", "accw_block", "clust"):
            a, b = m[m["class"] == "RA"][f"{k}|{s_}"].dropna(), m[m["class"] == "G4"][f"{k}|{s_}"].dropna()
            print(f"  {k:32s} {s_:6s} RA median {a.median():6.2f}  G4 median {b.median():6.2f}  MW p = {stats.mannwhitneyu(b, a).pvalue:.3g}")
    print(f"\n  regions overlapping a loci_focus DAR element: {(m.loci_focus_DARs != '').sum()} of {len(m)}")
    print(f"  wrote {O}/manual_phylop_scored.tsv, agreement.tsv")


if __name__ == "__main__":
    main()
