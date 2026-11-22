#!/usr/bin/env python
"""Does AR stand out among the predicted MAFB inputs for oscine site gain?

Every site in this analysis was called on lonStrDom2, an oscine. That alone
guarantees a positive mean gain for every matrix: a site has to be present in
the reference to be called, and the reference's close relatives inherit it.
The number that means something is therefore not AR's gain but AR's gain
RELATIVE TO the other predicted inputs, which carry the identical bias.

Reports, per matrix: mean gain with an element-clustered bootstrap CI, the
fraction of sites clearing the candidate threshold, and AR's rank. Then tests
AR against the pooled others (Mann-Whitney, and a permutation on element
labels so sites in one DAR are not treated as independent).
"""
import collections, json
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363"
                   "/phyloP_matched/loci_focus")
from span import F, DS, DE, elem_coords, in_domain, LABEL, LABEL_WIDE

M = f"{F}/maf_switch"
GAIN_MIN = 0.50
MN = {"full": "ARE-full (ref)", "jaspar": "JASPAR AR (ref)"}
REFS = set(MN.values())      # AR-pipeline matrices: shown for reference, not in the comparison set
ARNAMES = {"AR"}             # AR = the SCENIC+ matrix, scanned like every other regulator (= NR3C2)
FLAG = {"MAF", "STAT1", "LIN28B"}   # SCENIC+ motif evidence is another family's motif
META = FLAG
RGB = {"Glut-RA": "#e31a1c", "GABA-4": "#2c6e7f"}

el = elem_coords()

rec = []
for fn in ("data/site_gain.json", "data/tf_site_gain.json"):
    for g in json.load(open(f"{M}/{fn}")):
        c = el[g["elem"]]
        rec.append(dict(motif=MN.get(g["motif"], g["motif"]), cls=g["cls"],
                        elem=g["elem"], gain=g["f_osc"] - g["f_out"],
                        p_osc=g["p_osc"], n_osc=g["n_osc"],
                        p_out=g["p_out"], n_out=g["n_out"],
                        dom=in_domain(c)))

p = np.array([stats.fisher_exact([[r["p_osc"], r["n_osc"] - r["p_osc"]],
                                  [r["p_out"], r["n_out"] - r["p_out"]]])[1] for r in rec])
o = np.argsort(p)
q = np.empty(len(p))
q[o] = np.minimum(np.minimum.accumulate(
    (p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1], 1)
for r, qq in zip(rec, q):
    r["q"] = qq


def boot(rows, n=20000, seed=0):
    """Bootstrap the mean, resampling ELEMENTS so clustered sites stay together."""
    by = collections.defaultdict(list)
    for r in rows:
        by[r["elem"]].append(r["gain"])
    keys = list(by)
    if len(keys) < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    v = [np.array(by[k]) for k in keys]
    out = np.empty(n)
    for i in range(n):
        idx = rng.integers(0, len(v), len(v))
        out[i] = np.concatenate([v[j] for j in idx]).mean()
    return tuple(np.percentile(out, [2.5, 97.5]))


def table(rows, label):
    print(f"\n  === {label}  ({len(rows)} sites) ===")
    by = collections.defaultdict(list)
    for r in rows:
        by[r["motif"]].append(r)
    stat = []
    for m, v in by.items():
        g = np.array([x["gain"] for x in v])
        lo, hi = boot(v)
        stat.append((np.mean(g), m, len(v), np.median(g),
                     sum(1 for x in v if x["gain"] >= GAIN_MIN and x["q"] < 0.05),
                     len({x["elem"] for x in v}), lo, hi))
    stat.sort(reverse=True)
    print(f"    {'matrix':<16}{'n':>5}{'DARs':>6}{'mean':>8}{'95% CI':>18}"
          f"{'median':>8}{'cand':>6}")
    for mean, m, n, med, nc, nd, lo, hi in stat:
        star = " *" if m in ARNAMES else ("  m" if m in META else "")
        ci = f"[{lo:+.3f},{hi:+.3f}]" if np.isfinite(lo) else "        --       "
        print(f"    {m:<16}{n:>5}{nd:>6}{mean:>+8.3f}{ci:>18}{med:>+8.3f}{nc:>6}{star}")
    return stat


def contrast(rows, label):
    rows = [r for r in rows if r["motif"] not in REFS]   # references are not in the comparison
    ar = [r for r in rows if r["motif"] in ARNAMES]
    ot = [r for r in rows if not r["motif"] in ARNAMES]
    a = np.array([r["gain"] for r in ar]); b = np.array([r["gain"] for r in ot])
    if len(a) < 5 or len(b) < 5:
        print(f"    {label}: too few sites")
        return
    u = stats.mannwhitneyu(a, b, alternative="greater")
    # element-label permutation: shuffle the AR/other label over whole DARs
    by = collections.defaultdict(list)
    for r in rows:
        by[r["elem"]].append(r)
    keys = list(by)
    isar = {k: any(r["motif"] in ARNAMES for r in by[k]) for k in keys}
    nar = sum(isar.values())
    obs = a.mean() - b.mean()
    rng = np.random.default_rng(0); hit = 0; N = 20000
    for _ in range(N):
        pick = set(rng.choice(len(keys), nar, replace=False))
        x = np.concatenate([[r["gain"] for r in by[keys[i]]] for i in pick])
        y = np.concatenate([[r["gain"] for r in by[keys[i]]]
                            for i in range(len(keys)) if i not in pick])
        if x.mean() - y.mean() >= obs:
            hit += 1
    print(f"    {label}:  AR mean {a.mean():+.3f} (n={len(a)})   others {b.mean():+.3f}"
          f" (n={len(b)})   difference {obs:+.3f}")
    print(f"      Mann-Whitney (AR greater) p = {u.pvalue:.4f}"
          f"   element-permutation p = {(hit + 1) / (N + 1):.4f}")


def tiers(rows, label):
    """The phylogeny-aware version, on the avian363 tree: where does each site's
    principal oscine gain sit -- on the oscine stem, on the deep backbone just
    inside it (near-stem: nodes still holding >= 90% of oscines), shallower within
    oscines, or does the site predate oscines (inherited)?"""
    print(f"\n  === {label}: avian363 branch reconstruction ===")
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        if r.get("tier"): by[r["motif"]][r["tier"]] += 1
    T = ["stem", "near-stem", "within", "inherited", "none"]
    print(f"    {'matrix':<16}{'n':>5}{'stem':>6}{'near':>6}{'within':>8}{'inher':>7}{'none':>6}{'%stem+near':>12}")
    out = []
    key = lambda k: -((by[k]["stem"] + by[k]["near-stem"]) / max(1, sum(by[k].values())))
    for m in sorted(by, key=key):
        c = by[m]; n = sum(c.values())
        star = " *" if m in ARNAMES else ("  m" if m in META else "")
        print(f"    {m:<16}{n:>5}{c['stem']:>6}{c['near-stem']:>6}{c['within']:>8}{c['inherited']:>7}{c['none']:>6}"
              f"{100 * (c['stem'] + c['near-stem']) / n:>11.1f}%{star}")
        out.append(dict(motif=m, n=n, **{t: c[t] for t in T}))
    return out


tier_by = {}
for fn in ("trees363/counts_ar.json", "trees363/counts_tf.json"):
    for c in json.load(open(f"{M}/{fn}")):
        tier_by[(c["motif"], c["name"])] = c["tier"]
for r, g in zip(rec, [g for fn in ("data/site_gain.json", "data/tf_site_gain.json")
                      for g in json.load(open(f"{M}/{fn}"))]):
    r["tier"] = tier_by.get((g["motif"], g["name"]))

dom = table([r for r in rec if r["dom"]], f"DEFAULT SPAN -- {LABEL}")
win = table(rec, f"expanded span -- {LABEL_WIDE}")
td = tiers([r for r in rec if r["dom"]], f"DEFAULT SPAN -- {LABEL}")
tw = tiers(rec, f"expanded span -- {LABEL_WIDE}")
print("\n  AR versus the other predicted inputs")
contrast([r for r in rec if r["dom"]], "default span (domain)")
contrast(rec, "expanded span (+/-1.5 Mb)")

json.dump({"window": [{"motif": m, "n": n, "dars": nd, "mean": mean, "median": med,
                       "cand": nc, "lo": lo, "hi": hi}
                      for mean, m, n, med, nc, nd, lo, hi in win],
           "domain": [{"motif": m, "n": n, "dars": nd, "mean": mean, "median": med,
                       "cand": nc, "lo": lo, "hi": hi}
                      for mean, m, n, med, nc, nd, lo, hi in dom],
           "tier_window": tw, "tier_domain": td},
          open(f"{F}/results/compare_domain_mafb1p5mb/tf_gain_compare.json", "w"), indent=1)

# ---- figure -------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.2), sharex=True)
for ax, stat, ttl in ((axes[0], win, "MAFB $\\pm$1.5 Mb window"),
                      (axes[1], dom, "DHX35–TOP1 domain")):
    stat = sorted(stat)
    y = np.arange(len(stat))
    for i, (mean, m, n, med, nc, nd, lo, hi) in enumerate(stat):
        isar = m in ARNAMES
        col = "#e31a1c" if isar else ("#E59A9A" if m in REFS else ("#9aa6ad" if m in META else "#4c5a63"))
        if np.isfinite(lo):
            ax.plot([lo, hi], [i, i], color=col, lw=2.6 if isar else 1.6,
                    alpha=.95 if isar else .7, solid_capstyle="round", zorder=2)
        ax.plot([mean], [i], "o", ms=8 if isar else 5.5, color=col,
                mec="white", mew=1.1, zorder=3)
    ax.axvline(0, color="#b9c2c8", lw=1, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{m}  ({n})" for _, m, n, _, _, _, _, _ in stat], fontsize=9)
    for t, (_, m, _, _, _, _, _, _) in zip(ax.get_yticklabels(), stat):
        if m in ARNAMES:
            t.set_color("#e31a1c"); t.set_fontweight("bold")
        elif m in REFS:
            t.set_color("#E59A9A"); t.set_style("italic")
    ax.set_title(ttl, fontsize=11, pad=9)
    ax.set_xlabel("mean site gain  (oscine $-$ outgroup presence)", fontsize=9.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#e7ebee", lw=.8)
    ax.set_axisbelow(True)
fig.suptitle("Oscine site gain by matrix, predicted MAFB regulators\n"
             "all sites are reference-ascertained, so the baseline is positive for every matrix; "
             "red = AR (SCENIC+ matrix, = NR3C2); pink = AR reference matrices; grey = SCENIC+ evidence is another family's motif",
             fontsize=11.5, y=.985)
fig.tight_layout(rect=(0, 0, 1, .93))
for ext in ("pdf", "png"):
    fig.savefig(f"{F}/results/compare_domain_mafb1p5mb/fig_tf_gain.{ext}", dpi=190)
print(f"\n  wrote fig_tf_gain.pdf / .png and tf_gain_compare.json")
