#!/usr/bin/env python
"""Oscine site gain in the curated MAFB-domain lobes: do RA, G4 or shared DARs gain AR (or other regulator) sites?

Elements are the 110 curated lobes (results/mafb_domain/element_threshold_classes.tsv; RA, G4, shared = "both"),
each scanned over its scored extent on lonStrDom2. Matrices: the AR references (HOMER ARE-full = "full", JASPAR
MA0007.3 = "jaspar") and the SCENIC+ regulators (tf_scenic.meme; AR there is the NR3C1 matrix, so NR3C2 -- the
identical matrix -- is dropped). A site is a reference window scoring >= the matrix's 99.99th-percentile null
threshold (recount_cov.TH, the presence threshold), overlapping hits of one matrix merged to the best.

Presence per genome and ancestor is the corrected call (presence.py projection; missing below 80% motif coverage),
via tiers363.call on data/lobes_gap.maf (hal2maf over the lobes +/-100 bp, ancestors included). Per site:
  gain   oscine presence fraction - outgroup presence fraction (floors: >= 40 oscines, >= 60 outgroup covered)
  tier   avian363 placement of the principal oscine gain (tiers363.classify): stem / near-stem / within /
         inherited / none
Every site is ascertained on an oscine reference, so gain is positive for every matrix; compare classes and
matrices relative to one another, not to zero.

Tests, on all lobes and on lobes with osc and non-osc depth_min50 >= 0.8:
  mean gain by class, element-bootstrap CI; class pairs by a permutation of class labels over lobes (two-sided,
  Holm over the three pairs)
  stem + near-stem share by class (exact Fisher, Holm)
  AR minus other-regulator gain within each class (removes the shared ascertainment), compared between classes
  by the same lobe permutation

    cd maf_switch && python lobe_site_gain.py
Writes results/mafb_domain/element_motif_gain.tsv (one row per site) and element_motif_gain_summary.tsv.
"""
import collections, itertools, os, sys
import numpy as np, pandas as pd, pysam
from scipy import stats

sys.path.insert(0, "."); sys.path.insert(0, "..")
import presence as PR
from recount_cov import LO, RC, TH, CODE
from tiers363 import call, classify, OSC, OUT, STEM

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
FA = "/private/groups/colquittlab/assembly/lonStrDom2/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
MAF = "data/lobes_gap.maf"
FLOOR_OSC, FLOOR_OUT = 40, 60
DEPTH = 0.8
AR = {"full", "jaspar", "AR"}
DROP = {"NR3C2"}                       # same matrix as AR
FLAG = {"MAF", "STAT1", "LIN28B"}      # SCENIC+ motif is another family's
CLS = ["RA", "G4", "both"]
NAME = {"RA": "RA", "G4": "G4", "both": "shared"}
GAIN_MIN = 0.5
N = 20000
# MAFB eRegulon each SCENIC+ regulator belongs to (build_scenic_motifs.REG); the AR reference matrices are
# reported with RA but kept out of the eRegulon tests, so every TF enters once at one matrix
EREG = {"AR": "RA", "KLF12": "RA", "NFIA": "RA", "LIN28B": "RA", "STAT1": "RA", "FOSL2": "RA",
        "MAF": "G4", "SOX6": "G4", "TCF12": "G4", "ARX": "G4", "SATB1": "G4", "MEIS2": "G4"}
TIERS = ["stem", "near-stem", "within", "inherited", "none"]


def scan(seq, k):
    """-> [(start, score)] for windows clearing TH[k], overlapping hits merged to the best."""
    w = len(LO[k]); a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w: return []
    sl = np.lib.stride_tricks.sliding_window_view(a, w); ok = (sl < 4).all(1)
    idx = np.arange(w); sc = np.full(len(sl), -np.inf)
    s2 = np.where(sl < 4, sl, 0)
    sc[ok] = np.maximum(LO[k][idx, s2].sum(1), RC[k][idx, s2].sum(1))[ok]
    hits = [i for i in np.flatnonzero(sc >= TH[k])]
    out = []
    for i in hits:
        if out and i < out[-1][0] + w:
            if sc[i] > out[-1][1]: out[-1] = (i, sc[i])
        else:
            out.append((i, sc[i]))
    return out


def boot(d, col, seed=0):
    by = [g[col].values for _, g in d.groupby("element")]
    if len(by) < 3: return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    m = [np.concatenate([by[j] for j in rng.integers(0, len(by), len(by))]).mean() for _ in range(5000)]
    return tuple(np.percentile(m, [2.5, 97.5]))


def perm(d, stat, a, b, seed=0):
    """Two-sided permutation of class labels over the lobes of classes a and b."""
    d = d[d.cls.isin([a, b])]
    lobes = d.groupby("element").cls.first()
    obs = stat(d, set(lobes.index[lobes == a]), set(lobes.index[lobes == b]))
    if not np.isfinite(obs): return obs, np.nan
    rng = np.random.default_rng(seed); keys = lobes.index.values; na = (lobes == a).sum(); hit = 0
    for _ in range(N):
        p = rng.permutation(len(keys))
        v = stat(d, set(keys[p[:na]]), set(keys[p[na:]]))
        if abs(v) >= abs(obs) - 1e-12: hit += 1
    return obs, (hit + 1) / (N + 1)


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); out = np.empty(len(p)); run = 0
    for r, i in enumerate(o):
        run = max(run, min(1, (len(p) - r) * p[i])); out[i] = run
    return out


def mean_diff(d, A, B):
    x = d[d.element.isin(A)].gain; y = d[d.element.isin(B)].gain
    return x.mean() - y.mean() if len(x) and len(y) else np.nan


def ar_rel(d, A, B):
    """(AR - others) in A minus (AR - others) in B."""
    def rel(S):
        x = d[d.element.isin(S)]
        a, o = x[x.ar].gain, x[~x.ar].gain
        return a.mean() - o.mean() if len(a) and len(o) else np.nan
    return rel(A) - rel(B)


def build():
    el = pd.read_csv(f"{R}/element_threshold_classes.tsv", sep="\t")
    fa = pysam.FastaFile(FA)
    print("  loading alignment", flush=True)
    B = PR.load(MAF)
    mats = [k for k in LO if k not in DROP]
    rows = []
    for _, r in el.iterrows():
        seq = fa.fetch(r.chrom, r.ext_start, r.ext_end).upper()
        for k in mats:
            for i, sc in scan(seq, k):
                s = r.ext_start + i; e = s + len(LO[k])
                d = call(B, dict(s=s, e=e, motif=k))
                o = [g for g in OSC if d.get(g) is not None]; u = [g for g in OUT if d.get(g) is not None]
                rec = dict(element=r.element, cls=r.cls, thr_class=r.thr_class, chrom=r.chrom, start=s, end=e,
                           motif=k, score=round(float(sc), 3), ref_present=d.get("Lonchura_striata"),
                           n_osc=len(o), n_out=len(u), p_osc=sum(d[g] for g in o), p_out=sum(d[g] for g in u),
                           depth_ok=bool(r.osc_depth_min50 >= DEPTH and r.nonosc_depth_min50 >= DEPTH))
                c = classify(d)
                rec.update(tier=c["tier"], principal=c["principal"], gains=c["gains"], losses=c["losses"],
                           stem_state=d.get(STEM))
                rows.append(rec)
        print(f"    {r.element}: {sum(1 for x in rows if x['element'] == r.element)} sites", flush=True)
    d = pd.DataFrame(rows)
    d["floor_ok"] = (d.n_osc >= FLOOR_OSC) & (d.n_out >= FLOOR_OUT)
    d["f_osc"] = d.p_osc / d.n_osc.where(d.n_osc > 0)
    d["f_out"] = d.p_out / d.n_out.where(d.n_out > 0)
    d["gain"] = d.f_osc - d.f_out
    d["fisher_p"] = [stats.fisher_exact([[a, n - a], [b, m - b]], alternative="greater")[1] if ok else np.nan
                     for a, n, b, m, ok in zip(d.p_osc, d.n_osc, d.p_out, d.n_out, d.floor_ok)]
    d.to_csv(f"{R}/element_motif_gain.tsv", sep="\t", index=False, float_format="%.4g")
    return d


def report(d, label, summ):
    d = d[d.floor_ok].copy(); d["ar"] = d.motif.isin(AR)
    d["stem"] = d.tier.isin(["stem", "near-stem"])
    q = d.fisher_p.values; o = np.argsort(q); bh = np.empty(len(q))
    bh[o] = np.minimum(np.minimum.accumulate((q[o] * len(q) / np.arange(1, len(q) + 1))[::-1])[::-1], 1)
    d["cand"] = (d.gain >= GAIN_MIN) & (bh < 0.05)
    print(f"\n  ===== {label}: {d.element.nunique()} lobes with sites, {len(d)} sites =====")
    for grp, sub in (("AR (full + JASPAR + SCENIC AR)", d[d.ar]), ("other regulators", d[~d.ar]),
                     ("all matrices", d)):
        print(f"\n  -- {grp}")
        print(f"    {'class':<8}{'lobes':>6}{'sites':>7}{'mean gain':>11}{'95% CI':>18}{'stem+near':>12}{'cand':>7}")
        for c in CLS:
            x = sub[sub.cls == c]
            if not len(x): print(f"    {NAME[c]:<8}{0:>6}{0:>7}"); continue
            lo, hi = boot(x, "gain")
            print(f"    {NAME[c]:<8}{x.element.nunique():>6}{len(x):>7}{x.gain.mean():>+11.3f}"
                  f"{f'[{lo:+.3f},{hi:+.3f}]':>18}{f'{x.stem.sum()}/{len(x)}':>9}"
                  f"{100 * x.stem.mean():>3.0f}%{x.cand.sum():>4}")
            summ.append(dict(set=label, group=grp, cls=NAME[c], lobes=x.element.nunique(), sites=len(x),
                             mean_gain=x.gain.mean(), lo=lo, hi=hi, stem_near=int(x.stem.sum()), cand=int(x.cand.sum())))
        pairs = list(itertools.combinations(CLS, 2))
        pm = [perm(sub, mean_diff, a, b) for a, b in pairs]
        pf = [stats.fisher_exact([[sub[sub.cls == a].stem.sum(), (~sub[sub.cls == a].stem).sum()],
                                  [sub[sub.cls == b].stem.sum(), (~sub[sub.cls == b].stem).sum()]])[1]
              for a, b in pairs]
        hm, hf = holm([p for _, p in pm]), holm(pf)
        for (a, b), (obs, _), h1, h2 in zip(pairs, pm, hm, hf):
            print(f"    {NAME[a]} vs {NAME[b]}: gain diff {obs:+.3f}, lobe-permutation p = {h1:.3g};"
                  f"  stem+near Fisher p = {h2:.3g}  (Holm)")
    print("\n  -- AR relative to the other regulators in the same lobes' class (AR mean - others mean)")
    for c in CLS:
        x = d[d.cls == c]; a, o = x[x.ar].gain, x[~x.ar].gain
        print(f"    {NAME[c]:<8} AR {a.mean():+.3f} (n={len(a)})  others {o.mean():+.3f} (n={len(o)})"
              f"  AR - others {a.mean() - o.mean():+.3f}")
    pairs = list(itertools.combinations(CLS, 2)); pr = [perm(d, ar_rel, a, b) for a, b in pairs]
    for (a, b), (obs, _), h in zip(pairs, pr, holm([p for _, p in pr])):
        print(f"    {NAME[a]} vs {NAME[b]}: difference in (AR - others) {obs:+.3f}, lobe-permutation p = {h:.3g} (Holm)")
    print("\n  -- per matrix: mean gain (sites) by class")
    print(f"    {'matrix':<9}" + "".join(f"{NAME[c]:>15}" for c in CLS))
    for m, x in sorted(d.groupby("motif"), key=lambda t: -t[1].gain.mean()):
        cells = []
        for c in CLS:
            y = x[x.cls == c]; cells.append(f"{y.gain.mean():+.2f} ({len(y)})" if len(y) else "--")
        print(f"    {m + (' m' if m in FLAG else ''):<9}" + "".join(f"{v:>15}" for v in cells))


def ereg(d, label):
    """RA-eRegulon vs G4-eRegulon sites: what is gained, overall and in each DAR class.

    Main effect (RA vs G4 eRegulon TFs): exact permutation of the eRegulon labels over the 12 TFs (924 splits),
    the TF being the unit, since matrices differ in length and information and so in baseline gain.
    Interaction (does each eRegulon gain more in its own DAR class?): permutation of DAR class over lobes of
    the statistic (RA-reg - G4-reg gain in RA lobes) - (same in G4 lobes)."""
    d = d[d.floor_ok & d.motif.isin(EREG)].copy()
    d["ereg"] = d.motif.map(EREG); d["stem"] = d.tier.isin(["stem", "near-stem"])
    print(f"\n  ===== eRegulons, {label}: {len(d)} sites in {d.element.nunique()} lobes =====")
    print(f"    {'eRegulon':<9}{'DARs':<8}{'TFs':>4}{'lobes':>6}{'sites':>7}{'mean gain':>11}{'95% CI':>18}"
          + "".join(f"{t:>11}" for t in TIERS))
    for e in ("RA", "G4"):
        for c in CLS + ["all"]:
            x = d[d.ereg == e] if c == "all" else d[(d.ereg == e) & (d.cls == c)]
            lo, hi = boot(x, "gain")
            tc = x.tier.value_counts()
            print(f"    {e:<9}{NAME.get(c, c):<8}{x.motif.nunique():>4}{x.element.nunique():>6}{len(x):>7}"
                  f"{x.gain.mean():>+11.3f}{f'[{lo:+.3f},{hi:+.3f}]':>18}"
                  + "".join(f"{f'{tc.get(t, 0)} ({100 * tc.get(t, 0) / len(x):.0f}%)':>11}" for t in TIERS))
    # main effect, TF as the unit
    tf = d.groupby("motif").agg(gain=("gain", "mean"), stem=("stem", "mean"))
    tf["ereg"] = tf.index.map(EREG)
    for col, nm in (("gain", "mean gain"), ("stem", "stem+near share")):
        obs = tf[tf.ereg == "RA"][col].mean() - tf[tf.ereg == "G4"][col].mean()
        v = tf[col].values; null = [v[list(c)].mean() - np.delete(v, list(c)).mean()
                                    for c in itertools.combinations(range(len(v)), (tf.ereg == "RA").sum())]
        p = np.mean(np.abs(null) >= abs(obs) - 1e-12)
        print(f"    RA vs G4 eRegulon TFs ({nm}, per-TF means): diff {obs:+.3f}, exact TF-label permutation p = {p:.3g}")
    # site-level, Fisher on stem+near (sites not independent; descriptive)
    a, b = d[d.ereg == "RA"].stem, d[d.ereg == "G4"].stem
    print(f"    stem+near sites: RA eRegulons {a.sum()}/{len(a)}, G4 eRegulons {b.sum()}/{len(b)};"
          f" Fisher p = {stats.fisher_exact([[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]])[1]:.3g} (sites as independent)")
    # interaction
    def inter(x, A, B):
        def diff(S):
            y = x[x.element.isin(S)]; r, g = y[y.ereg == "RA"].gain, y[y.ereg == "G4"].gain
            return r.mean() - g.mean() if len(r) and len(g) else np.nan
        return diff(A) - diff(B)
    obs, p = perm(d, inter, "RA", "G4")
    print(f"    own-class interaction, gain: (RA-reg - G4-reg) in RA lobes minus in G4 lobes = {obs:+.3f},"
          f" lobe permutation p = {p:.3g}")
    d2 = d.assign(gain=d.stem.astype(float))
    obs, p = perm(d2, inter, "RA", "G4")
    print(f"    own-class interaction, stem+near share: {obs:+.3f}, lobe permutation p = {p:.3g}")


def main():
    d = build() if "--reuse" not in sys.argv else pd.read_csv(f"{R}/element_motif_gain.tsv", sep="\t")
    print(f"  {len(d)} sites; {d.floor_ok.sum()} clear the coverage floor; "
          f"{(d.ref_present != True).sum()} not called present in the reference projection")
    summ = []
    report(d, "all lobes", summ)
    report(d[d.depth_ok], "depth >= 0.8", summ)
    ereg(d, "all lobes"); ereg(d[d.depth_ok], "depth >= 0.8")
    pd.DataFrame(summ).to_csv(f"{R}/element_motif_gain_summary.tsv", sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
