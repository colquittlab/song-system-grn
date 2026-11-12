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


def centre(n, c, s, e, cls, SI, bws, size):
    """(centre, how): SnapATAC peak midpoint, else fragment summit within the element"""
    cl = UPCL.get(cls)
    if cl and c in SI[cl]:
        S, E = SI[cl][c]; i0 = max(0, np.searchsorted(S, s) - 50); i1 = np.searchsorted(S, e)
        ov = [(min(E[i], e) - max(S[i], s), i) for i in range(i0, i1) if E[i] > s and S[i] < e]
        if ov:
            i = max(ov)[1]; return (S[i] + E[i]) // 2, "snapatac"
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


def suffix(how):
    return ("_snapcenter" if how == "snapatac" else "") + ("_nosmooth" if K <= 1 else ("" if K == 11 else f"_smooth{K}"))


def run(span, how="element"):
    EL = elem_coords(); CLS = classes(); sfx = suffix(how)
    SI = snap_index() if how == "snapatac" else None
    where = span
    if span.startswith("set:"):                       # region BED from panelE_sets.py
        nm = span[4:]; where = "panelE"; sfx = f"_{nm}{sfx}"
        b = pd.read_csv(result(f"{nm}.bed", "panelE"), sep="\t", header=None, names=["c", "s", "e", "n", "k"])
        EL = {n: (c, s, e) for c, s, e, n in zip(b.c, b.s, b.e, b.n)}; CLS = dict(zip(b.n, b.k)); SI = None
    if span.startswith("set:"): names = list(EL)
    elif span == "domain": names = [n for n, c in EL.items() if in_domain(c)]
    elif span == "1p5mb": W = wide_elements(); names = [n for n in EL if n in W]
    else: names = list(EL)
    bws = {k: pyBigWig.open(p) for k, p in {**PHYLOP, **FRAG}.items()}
    size = {k: bw.chroms() for k, bw in bws.items()}
    base = {(p, t): [] for p in PAIRS for t in PHYLOP}       # base-level phyloP at high-diff bins
    rows = []
    for i, n in enumerate(sorted(names)):
        c, s, e = EL[n]
        if any(c not in size[k] for k in bws): continue
        mid, via = centre(n, c, s, e, CLS.get(n, ""), SI, bws, size) if SI else (s + (e - s) // 2, "element centre")
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
                         **boot(df, t)))
    S = pd.DataFrame(summ)
    S.to_csv(result(f"phylop_density_summary{sfx}.tsv", where), sep="\t", index=False, float_format="%.4g")
    print(S.to_string(index=False, float_format=lambda v: f"{v:.3g}"), flush=True)
    return df, base, S


def main():
    global K
    args = sys.argv[1:]; how = "element"
    if "--center" in args:
        i = args.index("--center"); how = args[i + 1]; del args[i:i + 2]
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
