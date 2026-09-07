"""Class-level summary of cross-species integration strength, finch x chicken vs finch x mouse.

Tests the write-up claim that all finch cell classes integrate well with chicken, while with
mouse non-neuronal cells integrate strongly, GABAergic neurons moderately and glutamatergic
neurons weakly. Reads the two full-suite hybrid-label composites (GSI + SAMap + CCA + SATURN):
  composite_scoring/results/gg_adult_hybrid/composite_calls.csv   (assemble_gg_adult_hybrid.py)
  composite_scoring/results/yao_adult_hybrid/composite_calls.csv  (assemble_yao_adult_hybrid.py)

Two per-cluster readouts, each summarised by finch cell class (Non-neuronal / GABAergic /
Glutamatergic, from the celltype_hybrid name prefix):
  top row     composite confidence -- the Zaremba magnitude channel evaluated at the rank-
              aggregated winner (composite_score.summarise); "is the best match strong?"
  bottom row  method agreement -- how many of the 4 methods independently pick the same
              reference label as their own top call; "do the methods agree on WHAT it is?"
Points are individual finch clusters; bars are class medians. Within-pair class differences
are tested by one-sided Mann-Whitney U (direction: Non-neuronal > GABA > Glut) and the
chicken-vs-mouse drop per class by paired Wilcoxon over shared finch clusters.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "font.size": 7, "axes.linewidth": 0.6})

BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring/results")
OUT = BASE / "class_summary"
OUT.mkdir(exist_ok=True)

PAIRS = [("gg_adult_hybrid", "Finch × chicken"), ("yao_adult_hybrid", "Finch × mouse")]
CLASSES = ["Non-neuronal", "GABAergic", "Glutamatergic"]
COLORS = {"Non-neuronal": "#1b9e77", "GABAergic": "#7570b3", "Glutamatergic": "#d95f02"}
INK, MUTED = "#222222", "#777777"


def cell_class(q):
    if q.startswith("GABA"):
        return "GABAergic"
    if q.startswith("Glut"):
        return "Glutamatergic"
    return "Non-neuronal"


def load(tag):
    D = pd.read_csv(BASE / tag / "composite_calls.csv", index_col=0)
    D["cls"] = [cell_class(q) for q in D.index]
    D["agree_frac"] = D.methods_agreeing_on_own_top / D.n_methods
    return D


def pstar(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


data = {tag: load(tag) for tag, _ in PAIRS}

# ---- per-class summary table -----------------------------------------------------------
rows = []
for tag, label in PAIRS:
    D = data[tag]
    for c in CLASSES:
        d = D[D.cls == c]
        rows.append(dict(pair=label, cell_class=c, n_clusters=len(d),
                         confidence_median=d.confidence.median(), confidence_mean=d.confidence.mean(),
                         frac_high_tier=(d.confidence_tier == "high").mean(),
                         frac_low_tier=(d.confidence_tier == "low").mean(),
                         agree_frac_mean=d.agree_frac.mean(),
                         frac_all_methods_agree=(d.methods_agreeing_on_own_top == d.n_methods).mean(),
                         frac_le1_method_agree=(d.methods_agreeing_on_own_top <= 1).mean()))
summ = pd.DataFrame(rows)
summ.to_csv(OUT / "class_summary_table.csv", index=False)

tests = []
for tag, label in PAIRS:
    D = data[tag]
    for met in ["confidence", "agree_frac"]:
        g = {c: D.loc[D.cls == c, met].values for c in CLASSES}
        for a, b in [("Non-neuronal", "GABAergic"), ("GABAergic", "Glutamatergic"), ("Non-neuronal", "Glutamatergic")]:
            p = mannwhitneyu(g[a], g[b], alternative="greater").pvalue
            tests.append(dict(pair=label, metric=met, test=f"{a} > {b}", p=p))
G, Y = data["gg_adult_hybrid"], data["yao_adult_hybrid"]
for c in CLASSES:
    idx = [q for q in G.index if G.cls[q] == c and q in Y.index]
    for met in ["confidence", "agree_frac"]:
        try:
            p = wilcoxon(G.loc[idx, met], Y.loc[idx, met], alternative="greater").pvalue
        except ValueError:
            p = np.nan
        tests.append(dict(pair="chicken vs mouse (paired)", metric=met, test=f"{c}: chicken > mouse", p=p))
tests = pd.DataFrame(tests)
tests.to_csv(OUT / "class_summary_tests.csv", index=False)
print(summ.round(3).to_string(index=False))
print()
print(tests.assign(p=tests.p.map(lambda x: f"{x:.2g}")).to_string(index=False))

# ---- figure -----------------------------------------------------------------------------
METRICS = [("confidence", "Composite confidence\n(match strength at winner)", (0, 1.0)),
           ("agree_frac", "Method agreement\n(fraction of 4 methods sharing top call)", (0, 1.0))]
fig, axes = plt.subplots(2, 2, figsize=(4.6, 4.2), sharey="row")
rng = np.random.default_rng(0)
for j, (tag, label) in enumerate(PAIRS):
    D = data[tag]
    for i, (met, ylabel, ylim) in enumerate(METRICS):
        ax = axes[i, j]
        for k, c in enumerate(CLASSES):
            v = D.loc[D.cls == c, met].values
            med = np.median(v) if met == "confidence" else np.mean(v)
            ax.bar(k, med, width=0.62, color=COLORS[c], alpha=0.28, linewidth=0, zorder=1)
            ax.hlines(med, k - 0.31, k + 0.31, color=COLORS[c], linewidth=1.6, zorder=3)
            if met == "agree_frac":   # discrete 0.25 steps: fan out horizontally, wrap rows of 6
                for val in np.unique(v):
                    n = (v == val).sum()
                    for r0 in range(0, n, 6):
                        m = min(6, n - r0)
                        xs = k + (np.arange(m) - (m - 1) / 2) * 0.075
                        ys = np.full(m, val) + 0.045 * (r0 // 6)
                        ax.scatter(xs, ys, s=9, color=COLORS[c], edgecolor="white",
                                   linewidth=0.4, zorder=4)
            else:
                xs = k + rng.uniform(-0.18, 0.18, len(v))
                ax.scatter(xs, v, s=9, color=COLORS[c], edgecolor="white", linewidth=0.4, zorder=4)
        # pairwise brackets: adjacent classes
        g = {c: D.loc[D.cls == c, met].values for c in CLASSES}
        for k, (a, b) in enumerate([("Non-neuronal", "GABAergic"), ("GABAergic", "Glutamatergic")]):
            p = mannwhitneyu(g[a], g[b], alternative="greater").pvalue
            y = 1.08
            ax.plot([k + 0.05, k + 0.05, k + 0.95, k + 0.95], [y - 0.02, y, y, y - 0.02],
                    color=INK, linewidth=0.6, clip_on=False)
            ax.text(k + 0.5, y + 0.01, pstar(p), ha="center", va="bottom", fontsize=6.5, color=INK)
        ax.set_ylim(0, 1.0)
        ax.set_xlim(-0.55, 2.55)
        ax.set_xticks(range(3))
        ns_ = [int((D.cls == c).sum()) for c in CLASSES]
        ax.set_xticklabels([f"Non-\nneuronal\n(n={ns_[0]})", f"GABA\n(n={ns_[1]})", f"Glut\n(n={ns_[2]})"],
                           fontsize=6.5)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=2, width=0.6, colors=INK, labelsize=6.5)
        ax.grid(axis="y", color="#e6e6e6", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        if met == "confidence":
            ax.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=2)
            if j == 1:
                ax.text(2.55, 0.5, " high-\n tier", va="center", ha="left", fontsize=5.5, color=MUTED)
        if j == 0:
            ax.set_ylabel(ylabel, fontsize=7)
        if i == 0:
            ax.set_title(label, fontsize=8, pad=16, color=INK)

fig.text(0.5, 0.01,
         "Points: finch clusters. Bars: class median (top) or mean (bottom).\n"
         "Brackets: one-sided Mann-Whitney U, left > right; * p<0.05, ** p<0.01, *** p<0.001, ns p>0.05.",
         ha="center", va="bottom", fontsize=5.5, color=MUTED)
fig.subplots_adjust(left=0.16, right=0.95, top=0.89, bottom=0.15, hspace=0.6, wspace=0.12)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"class_integration_summary.{ext}", dpi=300)
print(f"\nwrote {OUT}")
