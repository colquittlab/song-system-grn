#!/usr/bin/env python
"""Fig. 6E phyloP distribution, redone on the loci_focus DAR sets.

deeptools_mafb.qmd (comparative_genomics/deeptools) draws, for the MAFB-linked peaks, the
density of per-base phyloP at bases where cell-type accessibility differs strongly
(fragment-signal difference > 5), separately for Glut-RA − Glut-Arco-1 and
GABA-4-1 − GABA-1-1. Its peaks are SnapATAC2 cluster peaks filtered by ArchR p2g links to
MAFB, ArchR DARs (FDR < 0.1, Log2FC > 1) and 20–200 kb from the TSS. This script keeps the
measurement identical and swaps the peaks for the loci_focus DARs: the RA_over_Arco and
GABA4_over_GABA1 elements (dar_sets: FDR < 0.1, no fold-change cutoff; elements up in both
dropped) in a span.

Measurement, as in the notebook (computeMatrix reference-point, missing data as zero):
  ±1,000 bp around each element centre at 1 bp; phyloP and the two fragment differences
  each smoothed by an 11-bp centred rolling mean (edges left unsmoothed, as
  zoo::rollmean(fill = NA) then restored); keep bases with difference > 5.
Tracks: phyloP_avian363_3model_all ("avian all"), phyloP_mamAvi605_3model_all ("vert. all")
and phyloP_avian363_3model_birdAnc218 ("oscine (Anc218)", the core-oscine subtree, as in the
Fig. 6 phyloP tracks); fragment bigwigs Glut-RA, Glut-Arco-1, GABA-4-1, GABA-1-1 (bin 1, unsmoothed).

The Fig. 6E reading is a tail one: GABA-4 activity bases carry a heavier right tail (more
highly constrained bases) than Glut-RA activity bases. Reported per span and track:
  tail fraction  share of high-difference bases with smoothed phyloP > 2 and > 3, pooled over
                 bases (as the density pools them) -- descriptive
  per peak       each peak's tail fraction and mean phyloP over its high-difference bases;
                 RA vs G4 compared with Mann–Whitney (bases within a peak are not
                 independent, so the per-peak test is the one to quote)
  per base       the pooled (base-weighted) statistics -- share of high-difference bases
                 with phyloP > 1, 2, 3, 4 and the pooled mean -- G4 minus RA, with a 95%
                 interval and two-sided p from a peak-cluster bootstrap: peaks are drawn
                 with replacement and keep all their bases (in both contrasts), so the
                 within-element correlation of bases is respected. This is the test for
                 "more highly constrained bases overall", as the density shows it.

Centring. By default each window is centred on the 501-bp ArchR element. With
--center snapatac it is centred as Fig. 6E's are, on a SnapATAC2 cluster peak midpoint
(computeMatrix reference-point centre): an RA_over_Arco element takes the Glut-RA.bed peak
it overlaps most, a GABA4_over_GABA1 element the GABA-4-1.bed peak; an element with no
overlapping peak falls back to the summit of its up cluster's fragment signal (21-bp
rolling mean) within the element. Outputs then carry a _snapcenter suffix.
With --center summit each window is instead centred on the summit of the up cluster's fragment
signal (21-bp rolling mean) inside the SnapATAC2 peak the element overlaps most (or inside the
element if it overlaps none) -- the SnapATAC peak files carry no summit column, and their widths
(142-2,843 bp) put the midpoint well off the summit for wide peaks. Outputs carry _summit.

Track set. --tracks perbase swaps in the matched-neutral per-base phyloP tracks
(phyloP_matched/perbase/tracks; phyloP LRT CONACC against the matched neutral model):
"whole tree (matched)" and "oscine lineage, birdAnc222 (matched)". The latter is phyloP's
--subtree mode -- conservation or acceleration on the oscine clade relative to the rest of the
tree, i.e. oscine-lineage-specific constraint, not conservation within oscines. Outputs then
carry _perbase.

Window. --window W restricts each profile to a centred W-bp window (default 2000, i.e.
±1,000 bp as in the notebook); --window 500 is the element itself, --window 200 its central
±100 bp. Small windows also stop adjacent DARs' profiles overlapping. Outputs carry _w<W>.

Effect size. --lfc X keeps only elements whose own contrast's Log2FC exceeds X (Glut-RA vs
Glut-Arco-1 for RA_over_Arco, GABA-MGE-PVALB-2 vs GABA-LGE-1 for GABA4_over_GABA1; hybrid ArchR
tables exported by export_dar.R). G4 DARs are more strongly differential than RA DARs (domain
median Log2FC 2.9 vs 1.7), and a stronger difference puts more bases over the difference
threshold, so the cutoff also changes how the classes compare. Not applied to set: spans
(panelE_sets.py builds its own). Outputs carry _lfc<X>.

Matching. --match lfc pairs each RA element 1:1 with the G4 element nearest in own-contrast
Log2FC (greedy, without replacement, caliper 0.25), dropping the unpaired, so both classes
carry the same distribution of accessibility change. Outputs carry _matchlfc.

Tail mean. Every summary reports, for thresholds in TAILMEAN_Q, the mean phyloP of the
high-difference bases whose phyloP exceeds the threshold (pooled over bases), G4 minus RA with a
peak-cluster bootstrap interval -- how strongly constrained the constrained bases are,
alongside the tail fractions (how many there are).

Element means. Every summary also reports, for N in ELEM_N, each element's mean phyloP over its
high-difference bases, using only elements with >= N such bases; RA vs G4 compared per element
(difference of means with a 95% bootstrap interval over elements, and Mann-Whitney). Each
element counts once, unlike the base-weighted fractions.

Smoothing. --smooth K sets the rolling-mean window applied to phyloP and to the fragment
difference (default 11, as in the notebook); --smooth 1 uses raw per-base values for both, so
bases are also selected on unsmoothed accessibility. Outputs then carry _nosmooth.

Spans: domain (default), 1p5mb, gw (genome-wide reference), or set:<name> for a region BED
written by panelE_sets.py (results/panelE/<name>.bed: the rebuilt Fig. 6E set and the filter
ablation); set windows are centred on the region midpoint and outputs go to results/panelE/
with the set name in the file name.
    python phylop_density.py [--center snapatac] [--smooth K] [domain 1p5mb gw]
Writes results/<span>/phylop_density_{peaks,summary}.tsv and
results/compare_mafb_spans/fig_phylop_density.{png,pdf}.
"""
import sys, json
import numpy as np, pandas as pd, pyBigWig
from scipy import stats
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"], "pdf.fonttype": 42,
                            "font.size": 5, "axes.titlesize": 6, "axes.labelsize": 6,
                            "xtick.labelsize": 5, "ytick.labelsize": 5, "axes.linewidth": .4,
                            "xtick.major.width": .4, "ytick.major.width": .4})
import matplotlib.pyplot as plt

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0, F)
from span import elem_coords, in_domain, wide_elements, result

T = "/private/warm-archive/colquittlab/assembly/tracks_lark"
PHYLOP = {"avian all": f"{T}/phyloP/phyloP_avian363_3model_all.bw",
          "vert. all": f"{T}/phyloP/phyloP_mamAvi605_3model_all.bw",
          "oscine (Anc218)": f"{T}/phyloP/phyloP_avian363_3model_birdAnc218.bw"}
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
TRACK_SETS = {"3model": dict(PHYLOP),
              "perbase": {"whole tree (matched)": f"{PB}/wholetree_perbase_lsd2.bw",
                          "oscine lineage, birdAnc222 (matched)": f"{PB}/subtree222_perbase_lsd2.bw"}}
TRACKS = "3model"
LFC_MIN = None
MATCH = None
ELEM_N = (20, 50, 100)
TAILMEAN_Q = (0.0, 1.0, 2.0)
DART = "/private/warm-archive/colquittlab/assembly/tracks_lark/panelE_inputs"
FRAG = {k: f"{T}/fragment/{k}-fragment-bin1-smoothNone.bw" for k in ("Glut-RA", "Glut-Arco-1", "GABA-4-1", "GABA-1-1")}
PAIRS = {"Glut-RA − Glut-Arco-1": ("Glut-RA", "Glut-Arco-1"), "GABA-4-1 − GABA-1-1": ("GABA-4-1", "GABA-1-1")}
COL = {"Glut-RA − Glut-Arco-1": "#D6604D", "GABA-4-1 − GABA-1-1": "#4393C3"}
PAD, K, THR = 1000, 11, 5.0
TAILS = (2.0, 3.0)
BOOT_Q = (1.0, 2.0, 3.0, 4.0); NBOOT = 2000
SNAP = "/private/warm-archive/colquittlab/assembly/tracks_lark/snapatac_peaks"
UPCL = {"RA_over_Arco": "Glut-RA", "GABA4_over_GABA1": "GABA-4-1"}
SPAN_LABEL = {"domain": "DHX35–TOP1", "1p5mb": "MAFB ±1.5 Mb", "gw": "genome-wide"}


def rollmean(x, k=None):
    """centred rolling mean; the k//2 edge values keep their raw value (zoo fill = NA, restored)"""
    k = K if k is None else k
    if k <= 1: return x
    out = x.copy(); h = k // 2
    c = np.convolve(x, np.ones(k) / k, mode="valid")
    out[h:len(x) - h] = c
    return out


def vals(bw, c, s, e, size):
    s0, e0 = max(0, s), min(e, size)
    v = np.zeros(e - s, dtype=float)
    if e0 > s0:
        x = np.asarray(bw.values(c, s0, e0), dtype=float)
        v[s0 - s:e0 - s] = np.nan_to_num(x, nan=0.0)
    return v


def classes():
    cls = {}
    for l in open(f"{F}/allelem_cls.bed"):
        x = l.rstrip("\n").split("\t"); cls[x[3]] = x[4]
    return cls


def tail(df, bA, bB, t, q):
    """pooled and per-peak share of high-difference bases with phyloP > q, RA vs G4"""
    A, B = (df[f"tail{q:g}|{p}|{t}"].dropna().values for p in PAIRS)
    mw = stats.mannwhitneyu(B, A, alternative="two-sided") if len(A) and len(B) else None
    return {f"tail{q:g}_bases_RA": (bA > q).mean() if len(bA) else np.nan,
            f"tail{q:g}_bases_G4": (bB > q).mean() if len(bB) else np.nan,
            f"tail{q:g}_peakmean_RA": A.mean() if len(A) else np.nan,
            f"tail{q:g}_peakmean_G4": B.mean() if len(B) else np.nan,
            f"tail{q:g}_rank_biserial_G4_vs_RA": (2 * mw.statistic / (len(A) * len(B)) - 1) if mw else np.nan,
            f"tail{q:g}_p_peak": mw.pvalue if mw else np.nan}


def snap_index():
    """cluster -> chrom -> sorted (start, end) arrays of SnapATAC2 peaks"""
    out = {}
    for cl in set(UPCL.values()):
        d = pd.read_csv(f"{SNAP}/{cl}.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
        out[cl] = {c: (g.s.values, g.e.values) for c, g in d.sort_values(["c", "s"]).groupby("c")}
    return out


def centre(n, c, s, e, cls, SI, bws, size, summit=False):
    """(centre, how): SnapATAC peak midpoint (or, with summit, the up cluster's fragment summit
    inside that peak), else the fragment summit within the element"""
    cl = UPCL.get(cls)
    if cl and c in SI[cl]:
        S, E = SI[cl][c]; i0 = max(0, np.searchsorted(S, s) - 50); i1 = np.searchsorted(S, e)
        ov = [(min(E[i], e) - max(S[i], s), i) for i in range(i0, i1) if E[i] > s and S[i] < e]
        if ov:
            i = max(ov)[1]
            if not summit: return (S[i] + E[i]) // 2, "snapatac"
            f = rollmean(vals(bws[cl], c, int(S[i]), int(E[i]), size[cl][c]), 21)
            if f.max() > 0: return int(S[i]) + int(np.argmax(f)), "snapatac summit"
    if cl:
        f = rollmean(vals(bws[cl], c, s, e, size[cl][c]), 21)
        if f.max() > 0: return s + int(np.argmax(f)), "fragment summit"
    return s + (e - s) // 2, "element centre"


def boot(df, t, rng=None):
    """G4 - RA in base-weighted statistics, peak-cluster bootstrap (bases stay with their peak)"""
    rng = rng or np.random.default_rng(7)
    pa, pb = PAIRS                                   # RA pair, G4 pair
    nA, nB = df[f"nbins|{pa}"].values, df[f"nbins|{pb}"].values
    stat_cols = [("mean", f"sum|{pa}|{t}", f"sum|{pb}|{t}")] + \
                [(f"frac>{q:g}", f"n>{q:g}|{pa}|{t}", f"n>{q:g}|{pb}|{t}") for q in BOOT_Q]
    X = {k: (df[a].values, df[b].values) for k, a, b in stat_cols}
    obs = {k: xb.sum() / nB.sum() - xa.sum() / nA.sum() for k, (xa, xb) in X.items()}
    N = len(df); draws = {k: np.full(NBOOT, np.nan) for k in X}
    for i in range(NBOOT):
        ix = rng.integers(0, N, N); a, b = nA[ix].sum(), nB[ix].sum()
        if a == 0 or b == 0: continue                 # a draw with no bases in one contrast is undefined
        for k, (xa, xb) in X.items(): draws[k][i] = xb[ix].sum() / b - xa[ix].sum() / a
    out = {}
    for k in X:
        d = draws[k][np.isfinite(draws[k])]
        if len(d) < 100:
            lo = hi = p = np.nan
        else:
            lo, hi = np.percentile(d, [2.5, 97.5])
            p = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()) + 1 / len(d))
        out.update({f"base_{k}_RA": X[k][0].sum() / nA.sum(), f"base_{k}_G4": X[k][1].sum() / nB.sum(),
                    f"base_{k}_diff": obs[k], f"base_{k}_lo": lo, f"base_{k}_hi": hi, f"base_{k}_p": p})
    return out


def lfc_of():
    """element -> Log2FC in its own class's contrast"""
    out = {}; EL = elem_coords(); CLS = classes()
    for cl, con in (("RA_over_Arco", "ra"), ("GABA4_over_GABA1", "pvalb2_lge1")):
        d = pd.read_csv(f"{DART}/dar_{con}.tsv", sep="\t", usecols=["chrom", "end", "Log2FC"])
        m = dict(zip(zip(d.chrom, d.end), d.Log2FC))
        for n, (c, s, e) in EL.items():
            if CLS.get(n) == cl: out[n] = m.get((c, e), np.nan)
    return out


def match_lfc(names, L, CLS, cal=0.25):
    """1:1 nearest-Log2FC matching of RA to G4 elements, greedy without replacement"""
    ra = sorted((L[n], n) for n in names if CLS.get(n) == "RA_over_Arco" and np.isfinite(L.get(n, np.nan)))
    g4 = sorted((L[n], n) for n in names if CLS.get(n) == "GABA4_over_GABA1" and np.isfinite(L.get(n, np.nan)))
    gv = np.array([v for v, _ in g4]); used = np.zeros(len(g4), bool); keep = []
    # match the scarcer high-Log2FC RA elements first, so they are not crowded out
    for v, n in sorted(ra, reverse=True):
        if not len(gv): break
        j0 = np.searchsorted(gv, v); best = None
        for j in range(max(0, j0 - 200), min(len(gv), j0 + 200)):
            if used[j]: continue
            dlt = abs(gv[j] - v)
            if dlt <= cal and (best is None or dlt < best[0]): best = (dlt, j)
        if best: used[best[1]] = True; keep += [n, g4[best[1]][1]]
    return keep


def tail_means(df, t, rng=None):
    """G4 - RA in mean phyloP over high-difference bases with phyloP > q; peak-cluster bootstrap"""
    rng = rng or np.random.default_rng(13); pa, pb = PAIRS; out = {}; N = len(df)
    for q in TAILMEAN_Q:
        na, sa = df[f"tn>{q:g}|{pa}|{t}"].values, df[f"ts>{q:g}|{pa}|{t}"].values
        nb, sb = df[f"tn>{q:g}|{pb}|{t}"].values, df[f"ts>{q:g}|{pb}|{t}"].values
        mA = sa.sum() / na.sum() if na.sum() else np.nan; mB = sb.sum() / nb.sum() if nb.sum() else np.nan
        d = []
        for _ in range(NBOOT):
            ix = rng.integers(0, N, N); a, b = na[ix].sum(), nb[ix].sum()
            if a and b: d.append(sb[ix].sum() / b - sa[ix].sum() / a)
        d = np.array(d)
        ok = len(d) >= 100
        out.update({f"tailmean>{q:g}_RA": mA, f"tailmean>{q:g}_G4": mB, f"tailmean>{q:g}_nbases_RA": int(na.sum()),
                    f"tailmean>{q:g}_nbases_G4": int(nb.sum()), f"tailmean>{q:g}_diff": mB - mA,
                    f"tailmean>{q:g}_lo": np.percentile(d, 2.5) if ok else np.nan,
                    f"tailmean>{q:g}_hi": np.percentile(d, 97.5) if ok else np.nan,
                    f"tailmean>{q:g}_p": min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()) + 1 / len(d)) if ok else np.nan})
    return out


def elem_means(df, t):
    out = {}; rng = np.random.default_rng(11); pa, pb = PAIRS
    for N in ELEM_N:
        A = (df[f"sum|{pa}|{t}"] / df[f"nbins|{pa}"])[df[f"nbins|{pa}"] >= N].values
        B = (df[f"sum|{pb}|{t}"] / df[f"nbins|{pb}"])[df[f"nbins|{pb}"] >= N].values
        r = {f"elemN{N}_n_RA": len(A), f"elemN{N}_n_G4": len(B),
             f"elemN{N}_mean_RA": A.mean() if len(A) else np.nan, f"elemN{N}_mean_G4": B.mean() if len(B) else np.nan}
        if len(A) >= 3 and len(B) >= 3:
            dd = np.array([rng.choice(B, len(B)).mean() - rng.choice(A, len(A)).mean() for _ in range(NBOOT)])
            r.update({f"elemN{N}_diff": B.mean() - A.mean(), f"elemN{N}_lo": np.percentile(dd, 2.5),
                      f"elemN{N}_hi": np.percentile(dd, 97.5),
                      f"elemN{N}_p_mw": stats.mannwhitneyu(B, A, alternative="two-sided").pvalue})
        else:
            r.update({f"elemN{N}_diff": np.nan, f"elemN{N}_lo": np.nan, f"elemN{N}_hi": np.nan, f"elemN{N}_p_mw": np.nan})
        out.update(r)
    return out


def suffix(how):
    return ("_perbase" if TRACKS == "perbase" else "") + ("" if LFC_MIN is None else f"_lfc{LFC_MIN:g}") \
        + ("_matchlfc" if MATCH == "lfc" else "") + ("" if PAD == 1000 else f"_w{2 * PAD}") + ("_snapcenter" if how == "snapatac" else "_summit" if how == "summit" else "") + ("_nosmooth" if K <= 1 else ("" if K == 11 else f"_smooth{K}"))


def run(span, how="element"):
    EL = elem_coords(); CLS = classes(); sfx = suffix(how)
    SI = snap_index() if how in ("snapatac", "summit") else None
    where = span
    if span.startswith("set:"):                       # region BED from panelE_sets.py
        nm = span[4:]; where = "panelE"; sfx = f"_{nm}{sfx}"
        b = pd.read_csv(result(f"{nm}.bed", "panelE"), sep="\t", header=None, names=["c", "s", "e", "n", "k"])
        EL = {n: (c, s, e) for c, s, e, n in zip(b.c, b.s, b.e, b.n)}; CLS = dict(zip(b.n, b.k)); SI = None
    if span.startswith("set:"): names = list(EL)
    elif span == "domain": names = [n for n, c in EL.items() if in_domain(c)]
    elif span == "1p5mb": W = wide_elements(); names = [n for n in EL if n in W]
    else: names = list(EL)
    if (LFC_MIN is not None or MATCH) and not span.startswith("set:"):
        L = lfc_of()
        if LFC_MIN is not None:
            n0 = len(names); names = [n for n in names if L.get(n, np.nan) > LFC_MIN]
            print(f"  {span}: Log2FC > {LFC_MIN:g} keeps {len(names)} of {n0} elements", flush=True)
        if MATCH == "lfc":
            n0 = len(names); names = match_lfc(names, L, CLS)
            ra = [L[n] for n in names if CLS[n] == "RA_over_Arco"]; g4 = [L[n] for n in names if CLS[n] == "GABA4_over_GABA1"]
            print(f"  {span}: Log2FC matching keeps {len(ra)} RA + {len(g4)} G4 of {n0}; median Log2FC "
                  f"RA {np.median(ra):.2f}, G4 {np.median(g4):.2f}", flush=True)
    bws = {k: pyBigWig.open(p) for k, p in {**PHYLOP, **FRAG}.items()}
    size = {k: bw.chroms() for k, bw in bws.items()}
    base = {(p, t): [] for p in PAIRS for t in PHYLOP}       # base-level phyloP at high-diff bins
    rows = []
    for i, n in enumerate(sorted(names)):
        c, s, e = EL[n]
        if any(c not in size[k] for k in bws): continue
        mid, via = centre(n, c, s, e, CLS.get(n, ""), SI, bws, size, how == "summit") if SI else (s + (e - s) // 2, "element centre")
        a, b = mid - PAD, mid + PAD
        ph = {t: rollmean(vals(bws[t], c, a, b, size[t][c])) for t in PHYLOP}
        fr = {k: vals(bws[k], c, a, b, size[k][c]) for k in FRAG}
        r = dict(peak=n, chrom=c, start=s, end=e, cls=CLS.get(n, ""), centre=mid, centred_on=via,
                 shift=mid - (s + (e - s) // 2))
        for p, (x, y) in PAIRS.items():
            d = rollmean(fr[x] - fr[y]); m = d > THR
            r[f"nbins|{p}"] = int(m.sum())
            for t in PHYLOP:
                r[f"phylop|{p}|{t}"] = float(ph[t][m].mean()) if m.any() else np.nan
                for q in TAILS:
                    r[f"tail{q:g}|{p}|{t}"] = float((ph[t][m] > q).mean()) if m.any() else np.nan
                # per-peak sufficient statistics for the base-weighted bootstrap
                r[f"sum|{p}|{t}"] = float(ph[t][m].sum())
                for q in BOOT_Q: r[f"n>{q:g}|{p}|{t}"] = int((ph[t][m] > q).sum())
                for q in TAILMEAN_Q:
                    sel = ph[t][m][ph[t][m] > q]
                    r[f"tn>{q:g}|{p}|{t}"] = int(len(sel)); r[f"ts>{q:g}|{p}|{t}"] = float(sel.sum())
                if m.any(): base[(p, t)].append(ph[t][m])
        for t in PHYLOP: r[f"phylop_center100|{t}"] = float(ph[t][PAD - 50:PAD + 50].mean())
        rows.append(r)
        if span == "gw" and i % 5000 == 0: print(f"    {span}: {i}/{len(names)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(result(f"phylop_density_peaks{sfx}.tsv", where), sep="\t", index=False, float_format="%.4g")
    if SI is not None:
        print(f"  {span}: centred on " + ", ".join(f"{k} {v}" for k, v in df.centred_on.value_counts().items())
              + f"; |shift| median {df['shift'].abs().median():.0f} bp, 90th pct {df['shift'].abs().quantile(.9):.0f} bp", flush=True)
    summ = []
    for t in PHYLOP:
        A, B = (df[f"phylop|{p}|{t}"].dropna().values for p in PAIRS)
        mw = stats.mannwhitneyu(A, B, alternative="two-sided") if len(A) and len(B) else None
        rb = (2 * mw.statistic / (len(A) * len(B)) - 1) if mw else np.nan
        bA, bB = (np.concatenate(base[(p, t)]) if base[(p, t)] else np.array([]) for p in PAIRS)
        summ.append(dict(span=span.replace("set:", ""), centring=how, track=t, n_peaks=len(df),
                         n_peaks_RA=len(A), n_peaks_G4=len(B), n_bases_RA=len(bA), n_bases_G4=len(bB),
                         median_peakmean_RA=np.median(A) if len(A) else np.nan,
                         median_peakmean_G4=np.median(B) if len(B) else np.nan,
                         diff_medians=(np.median(A) - np.median(B)) if len(A) and len(B) else np.nan,
                         rank_biserial=rb, p_mannwhitney=mw.pvalue if mw else np.nan,
                         median_base_RA=np.median(bA) if len(bA) else np.nan,
                         median_base_G4=np.median(bB) if len(bB) else np.nan,
                         **{k: v for q in TAILS for k, v in tail(df, bA, bB, t, q).items()},
                         **boot(df, t), **tail_means(df, t), **elem_means(df, t)))
    S = pd.DataFrame(summ)
    S.to_csv(result(f"phylop_density_summary{sfx}.tsv", where), sep="\t", index=False, float_format="%.4g")
    print(S.to_string(index=False, float_format=lambda v: f"{v:.3g}"), flush=True)
    return df, base, S


def main():
    global K, TRACKS, PHYLOP, PAD, LFC_MIN, MATCH
    args = sys.argv[1:]; how = "element"
    if "--tracks" in args:
        i = args.index("--tracks"); TRACKS = args[i + 1]; del args[i:i + 2]
        PHYLOP = TRACK_SETS[TRACKS]
    if "--center" in args:
        i = args.index("--center"); how = args[i + 1]; del args[i:i + 2]
    if "--match" in args:
        i = args.index("--match"); MATCH = args[i + 1]; del args[i:i + 2]
    if "--lfc" in args:
        i = args.index("--lfc"); LFC_MIN = float(args[i + 1]); del args[i:i + 2]
    if "--window" in args:
        i = args.index("--window"); PAD = int(args[i + 1]) // 2; del args[i:i + 2]
    if "--smooth" in args:
        i = args.index("--smooth"); K = int(args[i + 1]); del args[i:i + 2]
    sfx = suffix(how)
    spans = args or ["domain", "1p5mb", "gw"]
    if all(sp.startswith("set:") for sp in spans): sfx = "_panelE" + sfx
    out = {sp: run(sp, how) for sp in spans}
    fig, axs = plt.subplots(len(spans), len(PHYLOP), figsize=(1.8 * len(PHYLOP) + .2, 1.05 * len(spans) + .3),
                            squeeze=False)
    for i, sp in enumerate(spans):
        df, base, S = out[sp]
        for j, t in enumerate(PHYLOP):
            ax = axs[i][j]
            for p in PAIRS:
                v = np.concatenate(base[(p, t)]) if base[(p, t)] else np.array([])
                if len(v) < 5: continue
                lo, hi = np.percentile(v, [.5, 99.5]); xs = np.linspace(lo, hi, 300)
                kd = stats.gaussian_kde(v)(xs)
                ax.plot(xs, kd, color=COL[p], lw=.7); ax.fill_between(xs, kd, color=COL[p], alpha=.15, lw=0)
            ax.axvline(0, color="#9AA3AE", lw=.3)
            ax.set_ylim(0, ax.get_ylim()[1] * 1.95)          # headroom for the annotation
            s = S[S.track == t].iloc[0]
            ax.set_title(f"{SPAN_LABEL.get(sp, sp.replace('set:', ''))} · {t}", loc="left", pad=2)
            ax.text(.98, .95, f"peaks {s.n_peaks_RA}/{s.n_peaks_G4}\nphyloP>2: {100*s['base_frac>2_RA']:.0f}% / {100*s['base_frac>2_G4']:.0f}%"
                    f"\nΔ {100*s['base_frac>2_diff']:+.1f} [{100*s['base_frac>2_lo']:+.1f}, {100*s['base_frac>2_hi']:+.1f}]"
                    + ("\nbootstrap p < 0.001" if s['base_frac>2_p'] < 1e-3 else f"\nbootstrap p = {s['base_frac>2_p']:.2g}"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=4.5, color="#4A525C")
            for k in ("top", "right"): ax.spines[k].set_visible(False)
            if i == len(spans) - 1: ax.set_xlabel("phyloP (per base)" if K <= 1 else f"phyloP ({K}-bp rolling mean)")
            if j == 0: ax.set_ylabel("density")
    h = [plt.Line2D([], [], color=c, lw=1, label=p) for p, c in COL.items()]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=5, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1), h_pad=.6, w_pad=.6)
    for ext in ("png", "pdf"):
        fig.savefig(result(f"fig_phylop_density{sfx}.{ext}", "compare_mafb_spans"), dpi=600)
    print(f"  wrote results/compare_mafb_spans/fig_phylop_density{sfx}.png/.pdf")


if __name__ == "__main__":
    main()
