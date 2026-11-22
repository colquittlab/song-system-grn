#!/usr/bin/env python
"""Does the visual rating read conservation relative to the local landscape?

A DAR on the shoulder of a strongly conserved neighbour can look weak to the eye even when its
absolute score is high. For each visually rated region, windows of the region's own width are
slid in 50-bp steps across +/-5 kb and scored like the region (whole tree: % of bases in 11-bp
blocks with mean phyloP > 1; oscine: % of bases with Anc222 phyloP > 1; both unweighted, so the
region and its neighbours are measured alike). Context scores:
  pctile   the region's rank among the windows within +/-5 kb
  fracmax  region score / best window within +/-2 kb (1 = the region is the local peak)
Each rating is then modelled on the absolute score (accessibility-weighted, best lobe for
multi-part regions; manual_phylop_lobes.py) plus a context score.

Writes results/manual_phylop/context.tsv.
"""
import numpy as np, pandas as pd, pyBigWig
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
O = f"{F}/results/manual_phylop"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
TR = {"whole": f"{PB}/wholetree_perbase_lsd2.bw", "oscine": f"{PB}/subtree222_perbase_lsd2.bw"}
FLANK, NEAR, STEP = 5000, 2000, 50
LAB = {0: "neutral", 1: "weak", 2: "moderate", 3: "strong"}


def score(v, t):
    return 100 * ((np.convolve(v, np.ones(11) / 11, "same") > 1) if t == "whole" else (v > 1)).mean()


def ols(y, X, names):
    X = np.column_stack([np.ones(len(y))] + X); b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b; s2 = r @ r / (len(y) - X.shape[1]); se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    p = 2 * stats.t.sf(abs(b / se), len(y) - X.shape[1]); r2 = 1 - (r @ r) / ((y - y.mean()) @ (y - y.mean()))
    return "  ".join(f"{n} {bb:+.2f} (p={pp:.2g})" for n, bb, pp in zip(["int"] + names, b, p)) + f"   R2 {r2:.2f}", r


def main():
    m = pd.read_csv(f"{O}/manual_phylop_scored.tsv", sep="\t")
    bws = {t: pyBigWig.open(p) for t, p in TR.items()}
    for t, bw in bws.items():
        pc, fm, rs = [], [], []
        for c, s, e in zip(m.chrom, m.start, m.end):
            s, e = int(s), int(e); w = e - s
            a0, b0 = max(0, s - FLANK), e + FLANK
            v = np.nan_to_num(np.array(bw.values(c, a0, b0), float))
            reg = score(v[s - a0:e - a0], t); rs.append(reg)
            win = [(x, score(v[x - a0:x - a0 + w], t)) for x in range(a0, b0 - w + 1, STEP)]
            sc = np.array([q for _, q in win]); near = np.array([q for x, q in win if abs(x - s) <= NEAR])
            pc.append(100 * (sc < reg).mean()); fm.append(reg / near.max() if near.max() > 0 else np.nan)
        m[f"ctx_{t}_region"] = rs; m[f"ctx_{t}_pctile"] = pc; m[f"ctx_{t}_fracmax"] = fm
    for t, vis in (("whole", "whole_r"), ("oscine", "oscine_r")):
        absolute = m[f"{t}_score_lobebest"]; z = lambda x: ((x - x.mean()) / x.std()).values
        print(f"\n== '{t}' rating")
        for k in ("pctile", "fracmax"):
            print(f"  Spearman with context {k:8s}: {stats.spearmanr(m[vis], m[f'ctx_{t}_{k}'])[0]:+.2f}")
        s0, r0 = ols(m[vis].values, [z(absolute)], ["absolute"])
        print(f"  absolute only:             {s0}")
        for k in ("pctile", "fracmax"):
            s1, r1 = ols(m[vis].values, [z(absolute), z(m[f'ctx_{t}_{k}'])], ["absolute", k])
            print(f"  absolute + {k:8s}:      {s1}")
        s2, _ = ols(m[vis].values, [z(absolute), z(m[f'ctx_{t}_fracmax']), (m["class"] == "G4").astype(float).values],
                    ["absolute", "fracmax", "G4"])
        print(f"  absolute + fracmax + G4:   {s2}")
        b = np.polyfit(absolute, m[vis], 1); res_abs = m[vis] - np.clip(np.polyval(b, absolute), 0, 3)
        X = np.column_stack([np.ones(len(m)), z(absolute), z(m[f'ctx_{t}_fracmax'])]); bb = np.linalg.lstsq(X, m[vis], rcond=None)[0]
        res_ctx = m[vis] - np.clip(X @ bb, 0, 3)
        m[f"{t}_resid_ctx"] = res_ctx
        print(f"  regions off by >= 1 level: absolute {(res_abs.abs() >= 1).sum()} -> with context {(res_ctx.abs() >= 1).sum()}")
        for reg in ("12751715-12752109", "12714605-12715179"):
            i = m.index[m.region.str.contains(reg)][0]
            print(f"    NC_042585.1:{reg} rated {LAB[int(m.loc[i, vis])]:8s} absolute {absolute[i]:5.1f}  local pctile {m.loc[i, f'ctx_{t}_pctile']:5.1f}  "
                  f"fraction of local max {m.loc[i, f'ctx_{t}_fracmax']:.2f}  residual {res_abs[i]:+.1f} -> {res_ctx[i]:+.1f}")
    cols = ["region", "class", "whole", "oscine"] + [c for c in m.columns if c.startswith("ctx_") or c.endswith("_resid_ctx")]
    m[cols].to_csv(f"{O}/context.tsv", sep="\t", index=False, float_format="%.3f")
    m.to_csv(f"{O}/manual_phylop_scored.tsv", sep="\t", index=False, float_format="%.3f")
    print(f"\n  wrote {O}/context.tsv")


if __name__ == "__main__":
    main()
