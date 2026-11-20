#!/usr/bin/env python
"""Report for lobe_gain_null.py: regulator site gain in the curated lobes against external expectations.

Measures, per site clearing the coverage floor:
  cand       gain >= 0.5 (oscine minus outgroup presence); at these sample sizes every such site has p < 1e-6
  stem+near  principal oscine gain on the oscine stem or the near-stem backbone (tiers363.classify)
Expectations:
  shuffled   the pooled rate of the matrix's own NSHUF column permutations (composition-matched)
  background the pooled rate of the random SCENIC+ matrices (motifs in general, same DNA)
Pooled and per-class comparisons: lobe bootstrap (lobes resampled with every set's sites in them, so the shared DNA
stays paired), 95% CI and a two-sided bootstrap p for regulator minus null. Per matrix: exact Poisson CI on the
observed count against each expected count. Branch rates (as tiers363.rate): gains per unit of covered branch
length on the stem + near-stem backbone, divided by the rate on the rest of the root-to-Lonchura path (fold = 1:
no excess at the origin of oscines), per set.
Writes results/mafb_domain/lobe_gain_null_{summary.tsv,per_matrix.tsv} and fig_lobe_gain_null.{pdf,png}.
"""
import collections, os, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                            "pdf.fonttype": 42, "ps.fonttype": 42})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
IN = os.environ.get("OUTF", f"{R}/lobe_gain_null_sites.tsv.gz")
REFS = {"full": "ARE-full (ref)", "jaspar": "JASPAR AR (ref)"}
EREG = {"AR": "RA", "KLF12": "RA", "NFIA": "RA", "LIN28B": "RA", "STAT1": "RA", "FOSL2": "RA",
        "MAF": "G4", "SOX6": "G4", "TCF12": "G4", "ARX": "G4", "SATB1": "G4", "MEIS2": "G4"}
FLAG = {"MAF", "STAT1", "LIN28B"}
CLS = ["RA", "G4", "both"]; NAME = {"RA": "RA", "G4": "G4", "both": "shared"}
NB = 2000
C_REG, C_AR, C_REF, C_SH, C_BG, INK, INK3, RULE = "#4C5A63", "#C8161A", "#E59A9A", "#7A9CC6", "#B7A27A", "#171B22", "#6B7380", "#DDE1E7"


def load():
    d = pd.read_csv(IN, sep="\t")
    d = d[d.floor_ok].copy()
    d["cand"] = d.gain >= 0.5
    d["stemnear"] = d.tier.isin(["stem", "near-stem"])
    d["grp"] = np.where(d.set == "regulator", np.where(d.matrix.isin(REFS), "ref", "regulator"),
                        np.where((d.set == "shuffled") & d.parent.isin(REFS), "ref_shuffled", d.set))
    return d


def boot_diff(d, col, a, b, mask=None, seed=0):
    """rate(a) - rate(b) and its lobe-bootstrap CI / two-sided p; a, b are grp names"""
    x = d if mask is None else d[mask]
    t = x.groupby(["element", "grp"])[col].agg(["sum", "count"]).unstack("grp").fillna(0)
    lob = t.index.values
    def rate(tt, g): return tt[("sum", g)].sum() / tt[("count", g)].sum() if tt[("count", g)].sum() else np.nan
    obs = (rate(t, a), rate(t, b))
    rng = np.random.default_rng(seed); diffs = []
    for _ in range(NB):
        tt = t.iloc[rng.integers(0, len(lob), len(lob))]
        diffs.append(rate(tt, a) - rate(tt, b))
    diffs = np.array(diffs)
    p = min(1, 2 * min((diffs <= 0).mean(), (diffs >= 0).mean()) + 1 / NB)
    lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
    return obs, obs[0] - obs[1], lo, hi, p, int(t[("count", a)].sum()), int(t[("count", b)].sum())


def poisson_ci(k):
    lo = stats.chi2.ppf(.025, 2 * k) / 2 if k > 0 else 0.0
    return lo, stats.chi2.ppf(.975, 2 * k + 2) / 2


def pois_p(k, mu):
    if mu <= 0: return np.nan
    return min(1, 2 * min(stats.poisson.cdf(k, mu), stats.poisson.sf(k - 1, mu)))


def branch(d, grp, mask=None, seed=0):
    x = d[d.grp == grp] if mask is None else d[(d.grp == grp) & mask]
    t = x.groupby("element")[["g_back", "bl_back", "g_rest", "bl_rest", "g_stem", "bl_stem"]].sum()
    def fold(tt, g, b): return (tt[g].sum() / tt[b].sum()) / (tt.g_rest.sum() / tt.bl_rest.sum())
    out = {}
    for lab, g, b in (("stem+near", "g_back", "bl_back"), ("stem", "g_stem", "bl_stem")):
        rng = np.random.default_rng(seed)
        bs = [fold(t.iloc[rng.integers(0, len(t), len(t))], g, b) for _ in range(NB)]
        out[lab] = (fold(t, g, b), *np.nanpercentile(bs, [2.5, 97.5]), int(t[g].sum()),
                    t[b].sum() * t.g_rest.sum() / t.bl_rest.sum())
    return out


def main():
    d = load()
    print(f"  sites clearing the floor: " + ", ".join(f"{k} {v}" for k, v in d.grp.value_counts().items())
          + f"; matrices: " + ", ".join(f"{k} {v}" for k, v in d.groupby('grp').matrix.nunique().items()))
    big = d[d.cand]
    print(f"  cand sites: max one-sided Fisher p {big.fisher_p.max():.2g}")
    summ, per = [], []
    for setlab, mask in (("all lobes", None), ("depth >= 0.8", d.depth_ok)):
        x = d if mask is None else d[mask]
        print(f"\n  ===== {setlab} =====")
        print(f"    {'':<20}{'sites':>7}{'mean gain':>11}{'cand':>9}{'stem+near':>11}")
        for g in ("regulator", "shuffled", "background", "ref", "ref_shuffled"):
            y = x[x.grp == g]
            print(f"    {g:<20}{len(y):>7}{y.gain.mean():>+11.3f}{100 * y.cand.mean():>8.1f}%{100 * y.stemnear.mean():>10.1f}%")
        for col, nm in (("cand", "cand"), ("stemnear", "stem+near"), ("gain", "mean gain")):
            for b in ("shuffled", "background"):
                (ra, rb), df, lo, hi, p, na, nb = boot_diff(d, col, "regulator", b, mask)
                print(f"    {nm:<10} regulators {ra:.3f} vs {b} {rb:.3f}: diff {df:+.3f} [{lo:+.3f},{hi:+.3f}], lobe-bootstrap p = {p:.3g}")
                summ.append(dict(set=setlab, measure=nm, cls="all", vs=b, reg=ra, null=rb, diff=df, lo=lo, hi=hi, p=p, n_reg=na, n_null=nb))
        print("    by DAR class (regulators vs background / shuffled):")
        for c in CLS:
            for col, nm in (("cand", "cand"), ("stemnear", "stem+near")):
                for b in ("background", "shuffled"):
                    m = (x.cls == c) if mask is None else (d.cls == c) & mask
                    (ra, rb), df, lo, hi, p, na, nb = boot_diff(d, col, "regulator", b, (d.cls == c) if mask is None else (d.cls == c) & mask)
                    summ.append(dict(set=setlab, measure=nm, cls=NAME[c], vs=b, reg=ra, null=rb, diff=df, lo=lo, hi=hi, p=p, n_reg=na, n_null=nb))
                    print(f"      {NAME[c]:<7}{nm:<10} reg {ra:.3f} (n={na}) vs {b} {rb:.3f} (n={nb}): {df:+.3f} [{lo:+.3f},{hi:+.3f}] p = {p:.3g}")
        print("    eRegulon x DAR class, cand / stem+near rate (background rate in the same class):")
        for c in CLS:
            y = x[x.cls == c]; bg = y[y.grp == "background"]
            cells = []
            for e in ("RA", "G4"):
                z = y[(y.grp == "regulator") & (y.matrix.map(EREG) == e)]
                cells.append(f"{e}-eReg {100 * z.cand.mean():.0f}% / {100 * z.stemnear.mean():.0f}% (n={len(z)})")
            print(f"      {NAME[c]:<7}" + "   ".join(cells) + f"   background {100 * bg.cand.mean():.0f}% / {100 * bg.stemnear.mean():.0f}%")
        print("    branch rates, stem + near-stem backbone (and stem alone) vs the rest of the path to Lonchura:")
        for g in ("regulator", "ref", "shuffled", "background"):
            br = branch(d, g, mask)
            a, s_ = br["stem+near"], br["stem"]
            print(f"      {g:<11} stem+near fold {a[0]:.2f} [{a[1]:.2f},{a[2]:.2f}] ({a[3]} gains, {a[4]:.1f} expected)"
                  f"   stem fold {s_[0]:.2f} [{s_[1]:.2f},{s_[2]:.2f}] ({s_[3]} gains, {s_[4]:.1f} expected)")
            for lab in ("stem+near", "stem"):
                summ.append(dict(set=setlab, measure=f"branch fold {lab}", cls="all", vs=g, reg=br[lab][0], lo=br[lab][1], hi=br[lab][2],
                                 n_reg=br[lab][3], null=br[lab][4]))
        # per matrix
        bgr = {col: x[x.grp == "background"][col].mean() for col in ("cand", "stemnear")}
        for m in sorted(x[x.set == "regulator"].matrix.unique()):
            y = x[x.matrix == m]; sh = x[(x.set == "shuffled") & (x.parent == m)]
            row = dict(set=setlab, matrix=REFS.get(m, m), n=len(y), lobes=y.element.nunique(), mean_gain=y.gain.mean(),
                       shuf_mean_gain=sh.gain.mean())
            for col in ("cand", "stemnear"):
                k = int(y[col].sum()); es = len(y) * sh[col].mean(); eb = len(y) * bgr[col]
                shr = sh.groupby("matrix")[col].mean()
                row.update({f"{col}_obs": k, f"{col}_exp_shuf": es, f"{col}_exp_bg": eb,
                            f"{col}_p_shuf": pois_p(k, es), f"{col}_p_bg": pois_p(k, eb),
                            f"{col}_frac_shuf_ge": (shr >= y[col].mean()).mean()})
            per.append(row)
        P = pd.DataFrame([r for r in per if r["set"] == setlab])
        print(f"    per matrix (obs / exp from own shuffles / exp from background; Poisson p):")
        for col, nm in (("cand", "cand"), ("stemnear", "stem+near")):
            print(f"      {nm}:")
            for r in P.sort_values(f"{col}_obs", ascending=False).itertuples(index=False):
                r = r._asdict()
                print(f"        {r['matrix']:<17}{r['n']:>4}  obs {r[col + '_obs']:>3}  shuf {r[col + '_exp_shuf']:>5.1f} (p {r[col + '_p_shuf']:.2g})"
                      f"  bg {r[col + '_exp_bg']:>5.1f} (p {r[col + '_p_bg']:.2g})  shuffles >= obs {100 * r[col + '_frac_shuf_ge']:.0f}%")
    pd.DataFrame(summ).to_csv(f"{R}/lobe_gain_null_summary.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(per).to_csv(f"{R}/lobe_gain_null_per_matrix.tsv", sep="\t", index=False, float_format="%.4g")
    figure(d, pd.DataFrame(per), pd.DataFrame(summ))


def figure(d, per, summ):
    fig = plt.figure(figsize=(12.6, 9.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.05], height_ratios=[1.6, 1], hspace=.38, wspace=.62,
                          left=.13, right=.98, top=.87, bottom=.07)
    P = per[(per.set == "all lobes") & (per.n >= 5)].copy()     # fewer than 5 sites: tables only (LIN28B, 1 site)
    P["is_ref"] = P.matrix.isin(REFS.values())
    pooled = {}
    x = d
    for col in ("cand", "stemnear"):
        r = x[x.grp == "regulator"]; sh = x[x.grp == "shuffled"]; bg = x[x.grp == "background"]
        pooled[col] = (int(r[col].sum()), len(r) * sh[col].mean(), len(r) * bg[col].mean(), len(r))
    for j, (col, ttl) in enumerate((("cand", "candidate gains  (gain ≥ 0.5)"), ("stemnear", "gains at or near the oscine stem"))):
        ax = fig.add_subplot(gs[0, j])
        Q = P.sort_values(["is_ref", f"{col}_obs"], ascending=[True, True])
        rows = [("all 12 regulators", *pooled[col], False)] + \
               [(r.matrix, r._asdict()[f"{col}_obs"], r._asdict()[f"{col}_exp_shuf"], r._asdict()[f"{col}_exp_bg"], r.n, r.is_ref)
                for r in Q.itertuples(index=False)]
        rows = rows[::-1]
        for i, (m, k, es, eb, n, isref) in enumerate(rows):
            lo, hi = poisson_ci(k)
            colr = C_AR if m == "AR" else (C_REF if isref else (INK if m.startswith("all") else C_REG))
            for off, e, mk, fc in ((.17, es, "o", colr), (-.17, eb, "D", "white")):
                if e <= 0: continue
                v = np.log2(max(k, .5) / e)
                ax.plot([np.log2(max(lo, .25) / e), np.log2(hi / e)], [i + off] * 2, color=colr, lw=1.4, alpha=.75, solid_capstyle="round")
                ax.plot([v], [i + off], mk, ms=5.2 if mk == "o" else 4.4, mfc=fc, mec=colr, mew=1.1, zorder=3)
        ax.axvline(0, color=INK3, lw=.9, zorder=0)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([f"{m}  {k}/{n}" for m, k, es, eb, n, _ in rows], fontsize=8.2)
        for t, (m, *_ , isref) in zip(ax.get_yticklabels(), rows):
            if m == "AR": t.set_color(C_AR); t.set_fontweight("bold")
            elif isref: t.set_color(C_REF); t.set_style("italic")
            elif m.startswith("all"): t.set_fontweight("bold")
        ax.set_xlabel("log$_2$ observed / expected", fontsize=9)
        ax.set_title(ttl, fontsize=10, loc="left")
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        ax.grid(axis="x", color=RULE, lw=.7); ax.set_axisbelow(True); ax.tick_params(labelsize=8)
    # panel: rates by DAR class
    ax = fig.add_subplot(gs[0, 2])
    G = [("regulator", "regulators", C_REG), ("shuffled", "shuffled", C_SH), ("background", "background", C_BG)]
    yl = []; y = 0
    for c in CLS + ["all"]:
        x = d if c == "all" else d[d.cls == c]
        for col, nm in (("cand", "cand"), ("stemnear", "stem+near")):
            for k, (g, lab, colr) in enumerate(G):
                z = x[x.grp == g]
                bs = []
                byl = z.groupby("element")[col].agg(["sum", "count"])
                rng = np.random.default_rng(k)
                for _ in range(1000):
                    t = byl.iloc[rng.integers(0, len(byl), len(byl))]; bs.append(t["sum"].sum() / t["count"].sum())
                lo, hi = np.percentile(bs, [2.5, 97.5]); v = z[col].mean()
                yy = y + (1 - k) * .24
                ax.plot([100 * lo, 100 * hi], [yy] * 2, color=colr, lw=1.6, solid_capstyle="round")
                ax.plot([100 * v], [yy], "o", ms=5, color=colr, mec="white", mew=.8, zorder=3)
            yl.append((y, f"{'all lobes' if c == 'all' else NAME[c]}  ·  {nm}")); y -= 1
        y -= .4
    ax.set_yticks([a for a, _ in yl]); ax.set_yticklabels([b for _, b in yl], fontsize=8.2)
    ax.set_xlabel("% of sites", fontsize=9); ax.set_title("rates by DAR class", fontsize=10, loc="left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis="x", color=RULE, lw=.7); ax.set_axisbelow(True); ax.tick_params(labelsize=8)
    ax.legend(handles=[Line2D([], [], color=c, marker="o", lw=1.6, label=l) for _, l, c in G], fontsize=7.8, frameon=False,
              loc="lower right")
    # panel: branch-rate folds
    ax = fig.add_subplot(gs[1, :2])
    bb = ax.get_position(); ax.set_position([bb.x0 + .08, bb.y0, bb.width - .08, bb.height])   # room for the row labels
    S = summ[(summ.measure.str.startswith("branch fold"))]
    labs = []; i = 0
    for setlab in ("all lobes", "depth >= 0.8"):
        for lab in ("stem+near", "stem"):
            for g, nm, colr in (("regulator", "regulators", C_REG), ("ref", "AR references", C_REF),
                                ("shuffled", "shuffled", C_SH), ("background", "background", C_BG)):
                r = S[(S.set == setlab) & (S.measure == f"branch fold {lab}") & (S.vs == g)]
                if r.empty: continue
                r = r.iloc[0]
                ax.plot([r.lo, r.hi], [i, i], color=colr, lw=1.6, solid_capstyle="round")
                ax.plot([r.reg], [i], "o", ms=5, color=colr, mec="white", mew=.8, zorder=3)
                ax.text(r.hi + .05 if np.isfinite(r.hi) else r.reg + .05, i, f"{int(r.n_reg)} gains / {r.null:.1f} exp", fontsize=7, va="center", color=INK3)
                labs.append((i, f"{'all' if setlab == 'all lobes' else 'depth ≥ 0.8'} · {lab} · {nm}")); i -= 1
            i -= .5
    ax.axvline(1, color=INK3, lw=.9)
    ax.set_yticks([a for a, _ in labs]); ax.set_yticklabels([b for _, b in labs], fontsize=7.6)
    ax.set_xlabel("gains per unit branch length at the oscine origin ÷ rate on the rest of the path to Lonchura", fontsize=9)
    ax.set_title("branch-rate test (avian363 tree)", fontsize=10, loc="left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis="x", color=RULE, lw=.7); ax.set_axisbelow(True); ax.tick_params(labelsize=8)
    # key
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.legend(handles=[Line2D([], [], color=C_REG, marker="o", lw=1.4, label="vs its own shuffles (composition-matched)"),
                       Line2D([], [], color=C_REG, marker="D", mfc="white", lw=1.4, label="vs random SCENIC+ matrices"),
                       Line2D([], [], color=C_AR, marker="o", lw=1.4, label="AR (SCENIC+ matrix)"),
                       Line2D([], [], color=C_REF, marker="o", lw=1.4, label="AR reference matrices (not pooled)")],
              loc="upper left", fontsize=8, frameon=False, title="per-matrix panels", title_fontsize=8.5)
    nb = d[d.grp == "background"].matrix.nunique(); ns = d[d.grp == "shuffled"].matrix.nunique()
    fig.suptitle("Oscine motif-site gain in the curated MAFB-domain lobes, against external expectations", fontsize=12.5,
                 fontweight="bold", y=.975)
    fig.text(.5, .942, f"110 lobes, sites clearing the coverage floor. Expectations from {ns} column-shuffled regulator matrices "
             f"and {nb} random SCENIC+ matrices, scanned and called identically.", ha="center", fontsize=8.6, color=INK3)
    fig.text(.5, .924, "Bars: exact Poisson 95% CI on the observed count (per matrix); lobe bootstrap 95% CI (rates, folds). "
             "Matrices with fewer than 5 sites are omitted from the per-matrix panels.", ha="center", fontsize=8.6, color=INK3)
    for ext in ("pdf", "png"):
        fig.savefig(f"{R}/fig_lobe_gain_null.{ext}", dpi=200)
    print(f"\n  wrote {R}/fig_lobe_gain_null.pdf / .png and lobe_gain_null_{{summary,per_matrix}}.tsv")


if __name__ == "__main__":
    main()
