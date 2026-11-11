#!/usr/bin/env python
"""Motif site gains per TF input, observed against expectation.

Every site in this analysis was called on lonStrDom2, an oscine, so every
matrix carries the same reference-ascertainment bias and every matrix is
expected to show gains. The expectation here is therefore empirical and
internal: the rate pooled across the regulator matrices (AR references excluded), applied to each
matrix's own number of tested sites.

    expected_m = n_m * (sum of gains over all matrices / sum of n)
    shown      = log2( observed_m / expected_m ), exact Poisson 95% CI on the
                 observed count

Two measures, because they ask different things:
  CANDIDATE GAINS   gain >= 0.50 and BH q < 0.05 in the span -- a frequency
                    asymmetry between oscines and the outgroup
  STEM/NEAR-STEM    the site's principal oscine gain, on the avian363 tree, sits on
                    the oscine stem or on the deep backbone just inside it (nodes
                    still holding >= 90% of oscines). Song likely arose through
                    several events near the base of oscines, so a gain a few nodes
                    inside the stem counts. Homoplasy elsewhere is not penalised.

Presence is the corrected call (maf_switch/presence.py): one reference-
projected sequence per genome; a genome with bases at fewer than 80% of the
motif's positions is missing, not absent.

Default span is the DHX35-TOP1 domain; the gene +/-1.5 Mb window is drawn
beside it as the expanded comparison.
"""
import collections, json, sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363"
                   "/phyloP_matched/loci_focus")
from span import F, elem_coords, in_domain, LABEL, LABEL_WIDE

M = f"{F}/maf_switch"
GAIN_MIN, Q_MAX = 0.50, 0.05
MN = {"full": "ARE-full (ref)", "jaspar": "JASPAR AR (ref)"}
REFS = set(MN.values())      # AR-pipeline matrices: shown for reference, not in the comparison set
ARNAMES = {"AR"}             # AR = the SCENIC+ matrix, scanned like every other regulator (= NR3C2)
# regulators whose SCENIC+ motif evidence is another family's motif
FLAG = {"MAF", "STAT1", "LIN28B"}
AR_C, META_C, TF_C, REF_C = "#C8161A", "#AEB6C2", "#4C5A63", "#E59A9A"
META = FLAG
INK, RULE = "#171B22", "#DDE1E7"

el = elem_coords()
tier_by = {}
for fn in ("trees363/counts_ar.json", "trees363/counts_tf.json"):
    for c in json.load(open(f"{M}/{fn}")):
        tier_by[(c["motif"], c["name"])] = c["tier"]

rec = []
for fn in ("data/site_gain.json", "data/tf_site_gain.json"):
    for g in json.load(open(f"{M}/{fn}")):
        rec.append(dict(motif=MN.get(g["motif"], g["motif"]),
                        dom=in_domain(el[g["elem"]]),
                        p_osc=g["p_osc"], n_osc=g["n_osc"],
                        p_out=g["p_out"], n_out=g["n_out"],
                        gain=g["f_osc"] - g["f_out"],
                        tier=tier_by.get((g["motif"], g["name"]))))


def bh(pv):
    pv = np.asarray(pv, float)
    o = np.argsort(pv)
    out = np.empty(len(pv))
    out[o] = np.minimum(np.minimum.accumulate(
        (pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 1)
    return out


def counts(rows, measure):
    """-> {matrix: (observed, n_sites)} with FDR computed within this span."""
    p = np.array([stats.fisher_exact([[r["p_osc"], r["n_osc"] - r["p_osc"]],
                                      [r["p_out"], r["n_out"] - r["p_out"]]])[1]
                  for r in rows])
    q = bh(p)
    out = collections.defaultdict(lambda: [0, 0])
    for r, qq in zip(rows, q):
        c = out[r["motif"]]
        c[1] += 1
        if measure == "cand":
            c[0] += int(r["gain"] >= GAIN_MIN and qq < Q_MAX)
        else:
            c[0] += int(r["tier"] in ("stem", "near-stem"))
    return {k: tuple(v) for k, v in out.items()}


def enrich(cnt):
    """log2 observed/expected under the pooled rate, with Poisson CI."""
    # the expectation pools the comparison set only; reference matrices are placed against it
    tot_k = sum(k for m, (k, _) in cnt.items() if m not in REFS)
    tot_n = sum(n for m, (_, n) in cnt.items() if m not in REFS)
    rate = tot_k / tot_n
    out = []
    for m, (k, n) in cnt.items():
        e = n * rate
        lo, hi = (stats.chi2.ppf(.025, 2 * k) / 2 if k else 0.0,
                  stats.chi2.ppf(.975, 2 * (k + 1)) / 2)
        out.append(dict(m=m, k=k, n=n, e=e,
                        v=np.log2((k + .5) / (e + .5)),
                        lo=np.log2((lo + .5) / (e + .5)),
                        hi=np.log2((hi + .5) / (e + .5))))
    out.sort(key=lambda d: d["v"])
    return out, rate


SPANS = (("domain", [r for r in rec if r["dom"]], f"{LABEL}  (default span)"),
         ("wide", rec, f"{LABEL_WIDE}  (expanded)"))
MEAS = (("cand", "candidate gains\ngain ≥ 0.50, q < 0.05"),
        ("stem", "gains at or near the oscine stem (avian363)\nstem, or a backbone node still holding \u226590% of oscines"))

fig, axes = plt.subplots(2, 2, figsize=(12.6, 10.4))
summary = {}
for col, (skey, rows, stitle) in enumerate(SPANS):
    for row, (mkey, mlab) in enumerate(MEAS):
        ax = axes[row][col]
        st, rate = enrich(counts(rows, mkey))
        summary[f"{skey}|{mkey}"] = st
        for i, d in enumerate(st):
            isar = d["m"] in ARNAMES; isref = d["m"] in REFS
            c = AR_C if isar else (REF_C if isref else (META_C if d["m"] in META else TF_C))
            ax.plot([d["lo"], d["hi"]], [i, i], color=c,
                    lw=2.8 if isar else 1.7, alpha=.95 if isar else .75,
                    solid_capstyle="round", zorder=2)
            ax.plot([d["v"]], [i], "o", ms=8.5 if isar else 6, color=c,
                    mec="white", mew=1.2, zorder=3)
            ax.text(ax.get_xlim()[1], i, "", fontsize=7)
        ax.axvline(0, color="#8C97A3", lw=1.1, zorder=1)
        ax.set_yticks(range(len(st)))
        ax.set_yticklabels([f'{d["m"]}   {d["k"]}/{d["n"]}' for d in st], fontsize=8.6)
        for t, d in zip(ax.get_yticklabels(), st):
            if d["m"] in ARNAMES:
                t.set_color(AR_C); t.set_fontweight("bold")
            elif d["m"] in REFS:
                t.set_color(REF_C); t.set_style("italic")
        ax.set_xlabel("log$_2$  observed / expected", fontsize=9.2)
        ax.tick_params(labelsize=8.4)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(RULE)
        ax.grid(axis="x", color="#EEF1F3", lw=.8)
        ax.set_axisbelow(True)
        if row == 0:
            ax.set_title(stitle, fontsize=11, color=INK, pad=24, fontweight="bold")
        ax.text(.015, 1.035, f"{mlab}    pooled rate {100 * rate:.0f}%",
                transform=ax.transAxes, fontsize=8.8, color="#5A6675", va="bottom")

fig.legend(handles=[Line2D([], [], marker="o", ls="none", color=AR_C, ms=8, label="AR (SCENIC+ matrix; = NR3C2)"),
                    Line2D([], [], marker="o", ls="none", color=REF_C, ms=6, label="AR reference matrix, not in the pool"),
                    Line2D([], [], marker="o", ls="none", color=TF_C, ms=6,
                           label="regulator matrix"),
                    Line2D([], [], marker="o", ls="none", color=META_C, ms=6,
                           label="SCENIC+ evidence is another family's motif")],
           loc="lower center", ncol=4, frameon=False, fontsize=8.6, bbox_to_anchor=(.5, .003))
fig.suptitle("Motif site gains per predicted MAFB input, against the pooled expectation",
             fontsize=13, fontweight="bold", color=INK, y=.983)
fig.text(.5, .952, "expectation is each matrix's own site count × the rate pooled over the regulator matrices "
                   "(AR references excluded); bars are exact Poisson 95% CIs. Labels give observed / sites tested.",
         ha="center", fontsize=9.2, color="#5A6675")
fig.tight_layout(rect=(0, .035, 1, .938))
for ext in ("pdf", "png"):
    fig.savefig(f"{F}/results/compare_domain_mafb1p5mb/fig_tf_gain_rank.{ext}", dpi=200)

json.dump(summary, open(f"{F}/results/compare_domain_mafb1p5mb/tf_gain_rank.json", "w"), indent=1, default=float)
for k, st in summary.items():
    print(f"\n  === {k} ===")
    print(f"    {'matrix':<16}{'obs':>5}{'sites':>7}{'exp':>8}{'log2 o/e':>10}")
    for d in sorted(st, key=lambda d: -d["v"]):
        star = " *" if d["m"] in ARNAMES else ""
        print(f"    {d['m']:<16}{d['k']:>5}{d['n']:>7}{d['e']:>8.1f}{d['v']:>+10.2f}{star}")
print("\n  wrote fig_tf_gain_rank.pdf / .png and tf_gain_rank.json")
