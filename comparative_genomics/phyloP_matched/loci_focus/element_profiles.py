#!/usr/bin/env python
"""Per-element constraint summaries and phyloP profiles for RA and G4 accessible elements.

Elements are SnapATAC2 cluster peaks (Glut-RA.bed for RA, GABA-4-1.bed for G4) in the span that
overlap an ArchR DAR of their own contrast (hybrid tables; FDR < 0.1, Log2FC > 0) -- the
"SnapATAC DARs" of the visual assessment. Unlike the loci_focus classes, elements up in both
contrasts are kept. Each peak is split into lobes at its up-cluster fragment summits: peaks
labelled bipartite / tripartite in the visual assessment get 2 / 3 lobes (the strongest summits
>= 150 bp apart); other peaks split where the 51-bp smoothed fragment signal has summits >= 150 bp
apart, each >= 30% of the strongest, with a dip between neighbours below 80% of the lower one
(at most three). Each lobe is centred on its own summit.

Per lobe, over the summit +/-100 bp (matched per-base phyloP tracks):
  whole_pct1 / osc_pct1     % bases with phyloP > 1           (whole tree / oscine lineage Anc222)
  whole_mean / osc_mean     mean phyloP
  whole_block / osc_block   % bases in 11-bp blocks with mean phyloP > 1 (continuity)
  whole_acc / osc_acc       % bases with phyloP < -1          (acceleration)
  log2fc                    own-contrast Log2FC of the overlapping DAR (hybrid ArchR tables)
and profiles over the summit +/-500 bp (11-bp rolling means): whole-tree phyloP, oscine phyloP,
and the up-cluster minus reference-cluster fragment difference.

    python element_profiles.py [domain|1p5mb]
Writes results/<span>/element_summaries.tsv and element_profiles_{whole,oscine,access}.tsv.gz.
"""
import sys
import numpy as np, pandas as pd, pyBigWig
from scipy.signal import find_peaks

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0, F)
from span import elem_coords, in_domain, wide_elements, result

T = "/private/warm-archive/colquittlab/assembly/tracks_lark"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
PHY = {"whole": f"{PB}/wholetree_perbase_lsd2.bw", "oscine": f"{PB}/subtree222_perbase_lsd2.bw"}
CLS = {"RA": ("RA_over_Arco", "Glut-RA", "Glut-Arco-1", "ra"), "G4": ("GABA4_over_GABA1", "GABA-4-1", "GABA-1-1", "pvalb2_lge1")}
CORE, PROF = 100, 500


def roll(x, k):
    return np.convolve(x, np.ones(k) / k, "same")


def vals(bw, c, s, e):
    size = bw.chroms()[c]; s0, e0 = max(0, s), min(e, size); v = np.zeros(e - s)
    if e0 > s0: v[s0 - s:e0 - s] = np.nan_to_num(np.array(bw.values(c, s0, e0), float))
    return v


def lobes(acc, forced=None):
    sm = roll(acc, 51)
    if sm.max() <= 0: return [int(len(acc) // 2)]
    if forced:                                       # labelled multi-part: strongest k summits >= 150 bp apart
        pk = []
        for i in np.argsort(-sm):
            if all(abs(i - j) >= 150 for j in pk): pk.append(int(i))
            if len(pk) == forced: break
        return sorted(pk)
    pk, _ = find_peaks(sm, distance=150, height=0.3 * sm.max())
    if len(pk) == 0: return [int(np.argmax(sm))]
    pk = sorted(pk, key=lambda i: -sm[i])[:3]; pk.sort(); keep = [pk[0]]
    for p in pk[1:]:
        dip = sm[keep[-1]:p + 1].min()
        if dip < 0.8 * min(sm[keep[-1]], sm[p]): keep.append(p)
        elif sm[p] > sm[keep[-1]]: keep[-1] = p
    return [int(i) for i in keep]


def main():
    span = sys.argv[1] if len(sys.argv) > 1 else "domain"
    lab = pd.read_csv(f"{F}/results/manual_phylop/manual_phylop.tsv", sep="\t", dtype=str).fillna("")
    FORCE = {}
    for reg, st in zip(lab.region, lab.structure.str.strip()):
        if st in ("bipartite", "tripartite"):
            a, b = reg.strip().split(":")[-1].replace(",", "").split("-")
            FORCE[(int(a) - 1, int(b))] = 2 if st == "bipartite" else 3
    EL = elem_coords(); cls = {}
    for l in open(f"{F}/allelem_cls.bed"):
        x = l.rstrip("\n").split("\t"); cls[x[3]] = x[4]
    W = wide_elements() if span == "1p5mb" else None
    keep = {n: c for n, c in EL.items() if (in_domain(c) if span == "domain" else n in W)}
    lfc = {}
    dar = {}
    for k, (cl, up, ref, con) in CLS.items():
        d = pd.read_csv(f"{T}/panelE_inputs/dar_{con}.tsv", sep="\t")
        d = d[(d.FDR < 0.1) & (d.Log2FC > 0)]; d["start"] = d.start - 1
        dar[k] = d[[(in_domain((c, s, e)) if span == "domain" else False) for c, s, e in zip(d.chrom, d.start, d.end)]]
    bw = {k: pyBigWig.open(p) for k, p in PHY.items()}
    fb = {x: pyBigWig.open(f"{T}/fragment/{x}-fragment-bin1-smoothNone.bw") for x in ("Glut-RA", "Glut-Arco-1", "GABA-4-1", "GABA-1-1")}
    rows, prof = [], {"whole": [], "oscine": [], "access": []}
    for k, (cl, up, ref, con) in CLS.items():
        sp = pd.read_csv(f"{T}/snapatac_peaks/{up}.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
        D = dar[k]
        for c, s, e in sp.itertuples(index=False):
            ov = D[(D.chrom == c) & (D.start < e) & (D.end > s)]
            if ov.empty: continue
            acc = vals(fb[up], c, s, e)
            ls = lobes(acc, FORCE.get((s, e)))
            for j, i in enumerate(ls):
                mid = s + i; a, b = mid - PROF, mid + PROF
                ph = {t: vals(bw[t], c, a, b) for t in PHY}
                df = roll(vals(fb[up], c, a, b) - vals(fb[ref], c, a, b), 11)
                core = slice(PROF - CORE, PROF + CORE)
                r = dict(element=f"{k}:{c}:{s + 1}-{e}" + (f"#{j + 1}" if len(ls) > 1 else ""), cls=k, chrom=c,
                         peak_start=s, peak_end=e, peak_width=e - s, n_lobes=len(ls), lobe=j + 1, summit=mid,
                         dars=",".join(f"{a}-{b}" for a, b in zip(ov.start, ov.end)), log2fc=ov.Log2FC.max(), fdr=ov.FDR.min(),
                         forced=(s, e) in FORCE)
                for t, lab in (("whole", "whole"), ("oscine", "osc")):
                    v = ph[t][core]
                    r[f"{lab}_pct1"] = 100 * (v > 1).mean(); r[f"{lab}_mean"] = v.mean()
                    r[f"{lab}_block"] = 100 * (roll(ph[t], 11)[core] > 1).mean(); r[f"{lab}_acc"] = 100 * (v < -1).mean()
                    prof[t].append(roll(ph[t], 11))
                prof["access"].append(df)
                rows.append(r)
    S = pd.DataFrame(rows)
    S.to_csv(result("element_summaries.tsv", span), sep="\t", index=False, float_format="%.3f")
    for t, M in prof.items():
        pd.DataFrame(np.array(M), index=S.element, columns=np.arange(-PROF, PROF)).to_csv(
            result(f"element_profiles_{t}.tsv.gz", span), sep="\t", float_format="%.3f")
    print(f"  {span}: {len(S)} lobes from {S.drop_duplicates(['cls', 'peak_start']).shape[0]} SnapATAC peaks; "
          + ", ".join(f"{k} {int((S.cls == k).sum())} lobes" for k in CLS)
          + f"; split peaks {S[S.n_lobes > 1].drop_duplicates(['cls', 'peak_start']).shape[0]}")
    got = {k: len(set(",".join(S[S.cls == k].dars).split(","))) for k in CLS}
    print("  ArchR DARs (FDR < 0.1, Log2FC > 0) covered by a SnapATAC peak: " + ", ".join(f"{k} {got[k]}/{len(dar[k])}" for k in CLS))


if __name__ == "__main__":
    main()
