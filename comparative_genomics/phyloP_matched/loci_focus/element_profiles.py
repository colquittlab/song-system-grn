#!/usr/bin/env python
"""Per-element constraint summaries and phyloP profiles for RA, G4 and shared accessible elements.

Elements are SnapATAC2 cluster peaks (Glut-RA.bed for RA, GABA-4-1.bed for G4) in the span that
overlap an ArchR DAR of their own contrast (hybrid tables; FDR < 0.1, Log2FC > 0) -- the
"SnapATAC DARs" of the visual assessment. Classes:
  RA, G4   peak of one class only
  both     an RA and a G4 peak overlapping each other (merged over their union), or a peak called
           in one cluster but overlapping DARs of both contrasts. For "both", accessibility is
           the RA and G4 up-cluster signals each scaled to its maximum and summed; Log2FC is the
           larger of the two; the accessibility profile is the mean of the two up-minus-reference
           differences.
Curation by eye, in files next to this script:
  excluded_peaks.bed    peaks dropped from the class named (not true or not differential peaks)
  manual_lobes.bed      lobes drawn on the browser; replace the automatic split, each centred on
                        its interval (the 1-bp pileup can be a plateau where the track shows lobes)
  class_override.bed    class set by eye
  lobe_count.bed        number of lobes by eye, split at the strongest summits >= 150 bp apart
Otherwise peaks are split into lobes at fragment summits: peaks labelled bipartite / tripartite in
the visual assessment get 2 / 3 lobes (strongest summits >= 150 bp apart); others split where
the 51-bp smoothed signal has summits >= 150 bp apart, each >= 30% of the strongest, with a dip
between neighbours below 80% of the lower one (at most three).

Per lobe, over the summit +/-100 bp (matched per-base phyloP tracks):
  whole_pct1 / osc_pct1     % bases with phyloP > 1           (whole tree / oscine lineage Anc222)
  whole_mean / osc_mean     mean phyloP
  whole_block / osc_block   % bases in 11-bp blocks with mean phyloP > 1 (continuity)
  whole_acc / osc_acc       % bases with phyloP < -1          (acceleration)
and over the lobe's full extent (ext_start-ext_end: the peak, cut between lobes at the smoothed
accessibility minimum, or the drawn interval), weighted by accessibility:
  whole_ext_pct1 / _block / _mean,  osc_ext_pct1 / _block / _mean
  acc_ra_ext / acc_g4_ext   mean Glut-RA / GABA-4-1 fragment signal over the extent (unweighted)
  osc_depth / nonosc_depth  alignment depth: mean fraction of oscine (145) / non-oscine (217) genomes
                            with a base aligned, over the extent (per-base tracks, domain_depth.py)
  *_depth_min50             the lowest 50-bp rolling mean of that depth within the extent
plus log2fc (own-contrast Log2FC of the overlapping DAR) and profiles over the summit +/-500 bp
(11-bp rolling means): whole-tree phyloP, oscine phyloP, accessibility difference.

    python element_profiles.py [domain|1p5mb]
Writes results/<span>/element_summaries.tsv and element_profiles_{whole,oscine,access}.tsv.gz.
"""
import sys
from pathlib import Path
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
    dar, dar_all = {}, {}
    for k, (cl, up, ref, con) in CLS.items():
        d = pd.read_csv(f"{T}/panelE_inputs/dar_{con}.tsv", sep="\t")
        d["start"] = d.start - 1; dar_all[k] = d
        d = d[(d.FDR < 0.1) & (d.Log2FC > 0)]
        dar[k] = d[[(in_domain((c, s, e)) if span == "domain" else False) for c, s, e in zip(d.chrom, d.start, d.end)]]
    bw = {k: pyBigWig.open(p) for k, p in PHY.items()}
    fb = {x: pyBigWig.open(f"{T}/fragment/{x}-fragment-bin1-smoothNone.bw") for x in ("Glut-RA", "Glut-Arco-1", "GABA-4-1", "GABA-1-1")}
    # candidate peaks per class: SnapATAC2 peaks of the up cluster overlapping an own-contrast DAR,
    # less peaks judged on inspection not to be true peaks (excluded_peaks.bed, next to this script)
    xb = [l.split("\t")[:4] for l in open(Path(__file__).with_name("excluded_peaks.bed")) if not l.startswith("#")]
    bad = lambda k, c, s, e: any(c == x and s < int(b) and int(a) < e and q in (k, "any") for x, a, b, q in xb)
    rd = lambda f: [l.rstrip("\n").split("\t") for l in open(Path(__file__).with_name(f)) if not l.startswith("#") and l.strip()]
    ml = [(x, int(a), int(b), q) for x, a, b, q, *_ in rd("manual_lobes.bed")]
    co = [(x, int(a), int(b), q) for x, a, b, q, *_ in rd("class_override.bed")]
    lc = [(x, int(a), int(b), int(q)) for x, a, b, q, *_ in rd("lobe_count.bed")]
    cand = {}
    for k, (cl, up, ref, con) in CLS.items():
        sp = pd.read_csv(f"{T}/snapatac_peaks/{up}.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
        D = dar[k]; cand[k] = []
        for c, s, e in sp.itertuples(index=False):
            ov = D[(D.chrom == c) & (D.start < e) & (D.end > s)]
            if not ov.empty and not bad(k, c, s, e): cand[k].append((c, s, e, ov))
    # RA and G4 peaks that overlap each other become one "both" element over their union
    pair = {(i, j) for i, a in enumerate(cand["RA"]) for j, b in enumerate(cand["G4"]) if a[0] == b[0] and a[1] < b[2] and b[1] < a[2]}
    pi, pj = {i for i, _ in pair}, {j for _, j in pair}
    assert len(pi) == len(pair) == len(pj), "RA/G4 peak overlaps are not one-to-one"
    # an unpaired peak that also overlaps a DAR of the other contrast (no SnapATAC peak called in that
    # cluster) is "both" too
    elems = []
    for k, o, used in (("RA", "G4", pi), ("G4", "RA", pj)):
        for i, (c, s, e, ov) in enumerate(cand[k]):
            if i in used: continue
            D = dar[o]; oo = D[(D.chrom == c) & (D.start < e) & (D.end > s)]
            q = [z for x, a, b, z in co if x == c and a < e and s < b]
            if q and q[0] == "both" and not len(oo):     # class set by eye: record the other contrast's DAR whatever its FDR
                A = dar_all[o]; oo = A[(A.chrom == c) & (A.start < e) & (A.end > s)]
            if q and q[0] != "both":
                if q[0] != k: A = dar_all[q[0]]; ov = A[(A.chrom == c) & (A.start < e) & (A.end > s)]
                elems.append((q[0], c, s, e, {q[0]: ov}))
            else: elems.append(("both", c, s, e, {k: ov, o: oo}, (s, e), (s, e)) if len(oo) else (k, c, s, e, {k: ov}))
    for i, j in sorted(pair):
        a, b = cand["RA"][i], cand["G4"][j]
        elems.append(("both", a[0], min(a[1], b[1]), max(a[2], b[2]), {"RA": a[3], "G4": b[3]}, (a[1], a[2]), (b[1], b[2])))
    rows, prof = [], {"whole": [], "oscine": [], "access": []}
    for el in elems:
        k, c, s, e, ovs = el[:5]
        ks = [x for x in ("RA", "G4") if x in ovs]
        if k == "both":
            # lobes on RA and G4 accessibility, each scaled to its own maximum, summed
            acc = sum(v / v.max() if v.max() > 0 else v for v in (vals(fb[CLS[x][1]], c, s, e) for x in ks))
            nf = max(FORCE.get(el[5], 0), FORCE.get(el[6], 0)) or None
        else:
            acc = vals(fb[CLS[k][1]], c, s, e); nf = FORCE.get((s, e))
        nf = next((n for x, a, b, n in lc if x == c and a < e and b > s), nf)
        man = sorted((max(a, s) - s, min(b, e) - s, q) for x, a, b, q in ml if x == c and a < e and b > s)
        if man:                                          # lobes drawn by eye: centred on each drawn interval
            ls = [(a + b) // 2 for a, b, q in man]; nf = len(man); ext = [(a, b) for a, b, q in man]
        else:
            ls = lobes(acc, nf); sm = roll(acc, 51)       # full extent: cut at the accessibility minimum between lobes
            cut = [p + int(np.argmin(sm[p:q])) for p, q in zip(ls[:-1], ls[1:])]
            ext = list(zip([0] + cut, cut + [e - s]))
        ph_ext = {t: vals(bw[t], c, s - 5, e + 5) for t in PHY}
        for j, i in enumerate(ls):
            mid = s + i; a, b = mid - PROF, mid + PROF
            ph = {t: vals(bw[t], c, a, b) for t in PHY}
            df = np.mean([roll(vals(fb[CLS[x][1]], c, a, b) - vals(fb[CLS[x][2]], c, a, b), 11) for x in ks], axis=0)
            core = slice(PROF - CORE, PROF + CORE)
            r = dict(element=f"{k}:{c}:{s + 1}-{e}" + (f"#{j + 1}" if len(ls) > 1 else ""), cls=k, chrom=c,
                     peak_start=s, peak_end=e, peak_width=e - s, n_lobes=len(ls), lobe=j + 1, summit=mid,
                     dars=";".join(f"{x}:" + ",".join(f"{p}-{q}" for p, q in zip(ovs[x].start, ovs[x].end)) for x in ks),
                     log2fc=max(ovs[x].Log2FC.max() for x in ks), fdr=min(ovs[x].FDR.min() for x in ks),
                     log2fc_ra=ovs["RA"].Log2FC.max() if "RA" in ovs else np.nan,
                     log2fc_g4=ovs["G4"].Log2FC.max() if "G4" in ovs else np.nan,
                     forced=nf is not None)
            for t, lab in (("whole", "whole"), ("oscine", "osc")):
                v = ph[t][core]
                r[f"{lab}_pct1"] = 100 * (v > 1).mean(); r[f"{lab}_mean"] = v.mean()
                r[f"{lab}_block"] = 100 * (roll(ph[t], 11)[core] > 1).mean(); r[f"{lab}_acc"] = 100 * (v < -1).mean()
                prof[t].append(roll(ph[t], 11))
                lo, hi = ext[j]; w = np.clip(acc[lo:hi], 0, None); w = w if w.sum() > 0 else np.ones(hi - lo)
                x = ph_ext[t][lo + 5:hi + 5]; blk = roll(ph_ext[t], 11)[lo + 5:hi + 5]
                r[f"{lab}_ext_pct1"] = 100 * (w * (x > 1)).sum() / w.sum(); r[f"{lab}_ext_block"] = 100 * (w * (blk > 1)).sum() / w.sum()
                r[f"{lab}_ext_mean"] = (w * x).sum() / w.sum()
            r["ext_start"], r["ext_end"] = s + ext[j][0], s + ext[j][1]
            for x, nm in (("Glut-RA", "acc_ra_ext"), ("GABA-4-1", "acc_g4_ext")):   # mean fragment signal over the extent
                r[nm] = vals(fb[x], c, s + ext[j][0], s + ext[j][1]).mean()
            prof["access"].append(df)
            rows.append(r)
    S = pd.DataFrame(rows)
    # alignment depth per lobe from the per-base clade depth tracks (domain_depth.py): mean over the
    # extent, and the lowest 50-bp rolling mean (a low-depth patch anywhere in the lobe)
    dbw = {x: pyBigWig.open(f"{F}/results/mafb_domain/depth_{x}.bw") for x in ("osc", "nonosc")}
    for x, bw_ in dbw.items():
        mean_, min_ = [], []
        for c, a, b in zip(S.chrom, S.ext_start, S.ext_end):
            v = np.array(bw_.values(c, int(a), int(b)), float)
            mean_.append(np.nanmean(v)); min_.append(np.nanmin(np.convolve(v, np.ones(50) / 50, "valid")) if len(v) >= 50 else np.nanmin(v))
        S[f"{x}_depth"], S[f"{x}_depth_min50"] = mean_, min_
    S.to_csv(result("element_summaries.tsv", span), sep="\t", index=False, float_format="%.3f")
    for t, M in prof.items():
        pd.DataFrame(np.array(M), index=S.element, columns=np.arange(-PROF, PROF)).to_csv(
            result(f"element_profiles_{t}.tsv.gz", span), sep="\t", float_format="%.3f")
    print(f"  {span}: {len(S)} lobes from {S.drop_duplicates(['cls', 'peak_start']).shape[0]} SnapATAC peaks; "
          + ", ".join(f"{k} {int((S.cls == k).sum())} lobes" for k in ("RA", "G4", "both"))
          + f"; split peaks {S[S.n_lobes > 1].drop_duplicates(['cls', 'peak_start']).shape[0]}")
    got = {k: len({d for x in S.dars for part in x.split(";") if part.startswith(k + ":") for d in part[len(k) + 1:].split(",")}) for k in CLS}
    print("  ArchR DARs (FDR < 0.1, Log2FC > 0) covered by a SnapATAC peak: " + ", ".join(f"{k} {got[k]}/{len(dar[k])}" for k in CLS))


if __name__ == "__main__":
    main()
