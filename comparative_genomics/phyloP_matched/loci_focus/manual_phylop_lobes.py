#!/usr/bin/env python
"""Score the multi-part (bipartite / tripartite) visually rated DARs lobe by lobe.

Lobes are defined from accessibility, not conservation, so the boundaries are not chosen by
the signal being scored: the up cluster's fragment signal (Glut-RA for RA, GABA-4-1 for G4;
51-bp rolling mean) is searched for as many peaks as the label says (2 or 3), at least 150 bp
apart, and the region is cut at the lowest point between neighbouring peaks. Each lobe is
scored with accessibility weighting, as for the single-part regions (manual_phylop_compare.py):
  whole   accessibility-weighted % of bases in 11-bp blocks with mean phyloP > 1 (whole tree)
  oscine  accessibility-weighted % of bases with oscine-lineage (Anc222) phyloP > 1
and each region is summarised by its best lobe and its mean over lobes. Agreement with the
visual ratings is then recomputed with the lobe scores in place of the full-region scores.

Writes results/manual_phylop/{lobes.tsv, lobes.bed}.
"""
import numpy as np, pandas as pd, pyBigWig
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
O = f"{F}/results/manual_phylop"
FR = "/private/warm-archive/colquittlab/assembly/tracks_lark/fragment"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
UP = {"RA": "Glut-RA", "G4": "GABA-4-1"}
W, OC = "whole tree (matched)", "oscine lineage Anc222 (matched)"
NLOBE = {"bipartite": 2, "tripartite": 3}
LAB = {0: "neutral", 1: "weak", 2: "moderate", 3: "strong"}


def vals(bw, c, s, e):
    return np.nan_to_num(np.array(bw.values(c, s, e), float))


def lobes(acc, k, gap=150):
    """cut points splitting acc into k lobes at the minima between its k strongest peaks"""
    sm = np.convolve(acc, np.ones(51) / 51, "same"); order = np.argsort(-sm); peaks = []
    for i in order:
        if all(abs(i - j) >= gap for j in peaks): peaks.append(int(i))
        if len(peaks) == k: break
    peaks.sort()
    cuts = [a + int(np.argmin(sm[a:b])) for a, b in zip(peaks[:-1], peaks[1:])]
    return [0] + cuts + [len(acc)], peaks


def accw(v, w, block=False):
    x = (np.convolve(v, np.ones(11) / 11, "same") > 1) if block else (v > 1)
    return 100 * (x * w).sum() / w.sum() if w.sum() else np.nan


def main():
    m = pd.read_csv(f"{O}/manual_phylop_scored.tsv", sep="\t").fillna({"structure": ""})
    lab = pd.read_csv(f"{O}/manual_phylop.tsv", sep="\t", dtype=str).fillna("")
    m["structure"] = lab.structure.str.strip().values            # labels as now in the sheet
    bw = {"whole": pyBigWig.open(f"{PB}/wholetree_perbase_lsd2.bw"), "oscine": pyBigWig.open(f"{PB}/subtree222_perbase_lsd2.bw")}
    fb = {k: pyBigWig.open(f"{FR}/{v}-fragment-bin1-smoothNone.bw") for k, v in UP.items()}
    rows, bed = [], []
    for i, r in m[m.structure.isin(NLOBE)].iterrows():
        c, s, e = r.chrom, int(r.start), int(r.end)
        acc = vals(fb[r["class"]], c, s, e); vw, vo = vals(bw["whole"], c, s, e), vals(bw["oscine"], c, s, e)
        cuts, peaks = lobes(acc, NLOBE[r.structure])
        sc = []
        for j, (a, b) in enumerate(zip(cuts[:-1], cuts[1:])):
            w = acc[a:b]
            sc.append(dict(lobe=j + 1, start=s + a, end=s + b, whole=accw(vw[a:b], w, block=True), oscine=accw(vo[a:b], w)))
            bed.append((c, s + a, s + b, f"{r.region.strip()}|lobe{j + 1}|whole {sc[-1]['whole']:.0f}|oscine {sc[-1]['oscine']:.1f}",
                         0, ".", s + a, s + b, "199,82,27" if j % 2 == 0 else "44,110,127"))
        for t in ("whole", "oscine"):
            m.loc[i, f"lobe_best_{t}"] = max(x[t] for x in sc); m.loc[i, f"lobe_mean_{t}"] = np.mean([x[t] for x in sc])
        rows.append(dict(region=r.region.strip(), cls=r["class"], structure=r.structure,
                         rated_whole=LAB[int(r.whole_r)], rated_oscine=LAB[int(r.oscine_r)],
                         full_whole=r[f"{W}|accw_block"], full_oscine=r[f"{OC}|accw_pct1"],
                         **{f"lobe{x['lobe']}": f"chr20:{x['start'] + 1:,}-{x['end']:,} (whole {x['whole']:.1f}, oscine {x['oscine']:.1f})" for x in sc}))
    L = pd.DataFrame(rows); L.to_csv(f"{O}/lobes.tsv", sep="\t", index=False, float_format="%.1f")
    pd.DataFrame(bed).to_csv(f"{O}/lobes.bed", sep="\t", header=False, index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    print(L.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    # agreement with lobe scores substituted for the multi-part regions
    for t, vis, full in (("whole", "whole_r", f"{W}|accw_block"), ("oscine", "oscine_r", f"{OC}|accw_pct1")):
        out = [f"\n'{t}': Spearman rho  full region {stats.spearmanr(m[vis], m[full])[0]:+.3f}"]
        for how in ("best", "mean"):
            x = m[full].where(m[f"lobe_{how}_{t}"].isna(), m[f"lobe_{how}_{t}"])
            m[f"{t}_score_lobe{how}"] = x
            out.append(f"{how} lobe {stats.spearmanr(m[vis], x)[0]:+.3f}")
        print("   ".join(out))
        for how in ("best", "mean"):
            x = m[f"{t}_score_lobe{how}"]; b = np.polyfit(x, m[vis], 1); res = m[vis] - np.clip(np.polyval(b, x), 0, 3)
            mp = m.structure.isin(NLOBE)
            print(f"   {how} lobe: multi-part residuals " + ", ".join(f"{rr.region.strip().split(':')[-1]} {q:+.1f}" for (_, rr), q in zip(m[mp].iterrows(), res[mp]))
                  + f"   | regions off by >= 1 level: {(res.abs() >= 1).sum()}")
    m.to_csv(f"{O}/manual_phylop_scored.tsv", sep="\t", index=False, float_format="%.3f")


if __name__ == "__main__":
    main()
