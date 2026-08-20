"""Cluster x k.anchor heatmap of the reciprocity-weighted RPCA sweep (hybrid-labeled
glut-only finch x excitatory chicken comparison) -- companion to
plot_sweep_percluster_glutonly_hybrid.py's titration-curve view. Where that plot shows
group means with individual lines as low-alpha context, this shows every cluster's own
row explicitly.

Rows use a FIXED biological order (per explicit user request), not hierarchical
clustering by response profile: the five song nuclei in HVCra, HVCx, LMANco, LMANsh, RA
order, then the non-song DACH2 clusters in numeric order, then the non-song CACNA1H
clusters in numeric order.

Reuses the *_reciprocal.csv outputs already written by plot_sweep_percluster_glutonly_hybrid.py
rather than recomputing -- those already exclude Glut-DACH2-HVCra-Pre per the user's
earlier request, and already classify Glut-DACH2-LMANsh as song.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from position_palette import POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
KA = [2, 5, 10, 20, 30, 50, 75, 100]
OUT = Path("/private/groups/colquittlab/saturn/zaremba_composite/results/composite_gg_glutonly_hybrid")

VARIANTS = [
    ("rpca_kfNA", "RPCA (hybrid labels), log-norm, k.filter=NA"),
    ("rpca_kf200", "RPCA (hybrid labels), log-norm, k.filter=200"),
]

# Fixed biological row order (user-specified): song nuclei in HVCra, HVCx, LMANco,
# LMANsh, RA order, then non-song DACH2 clusters and non-song CACNA1H clusters each in
# numeric order.
ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]

loaded = {}
for tag, title in VARIANTS:
    D = pd.read_csv(OUT / f"sweep_percluster_{tag}_reciprocal.csv", index_col=0)
    group = D["group"]
    M = D[[str(k) for k in KA]] if str(KA[0]) in D.columns else D[[c for c in D.columns if c != "group"]]
    M.columns = [int(c) for c in M.columns]
    M = M[KA]
    loaded[tag] = (M, group)

# Fixed vmax shared across BOTH k.filter panels (project convention, see
# plot_rank_heatmap.py): a paler panel must mean weaker integration, not just a smaller
# max within that panel alone.
vmax = max(M.values.max() for M, _ in loaded.values())

# Row pitch derived directly from the 6pt label size (not a fixed figure height
# regardless of row count) so cell height tracks the label text tightly instead of
# leaving tall blank-looking rows around each label.
LABEL_PT = 6
n_rows = len(ROW_ORDER)
row_pitch_in = LABEL_PT / 72.0 * 1.15
core_h = row_pitch_in * n_rows
fig_h = core_h + 1.55       # + bottom (x ticks/xlabel) + per-axis title + suptitle/legend margin
n_cols = len(KA)
core_w = row_pitch_in * n_cols     # same pitch as rows -> square cells
fig_w = 2 * (core_w + 1.9) + 0.6   # + per-panel margin (y labels + colorbar) + inter-panel gap

fig, axes = plt.subplots(1, len(VARIANTS), figsize=(fig_w, fig_h))
if len(VARIANTS) == 1:
    axes = [axes]

for ax, (tag, title) in zip(axes, VARIANTS):
    M, group = loaded[tag]

    order = [c for c in ROW_ORDER if c in M.index]
    missing = set(M.index) - set(order)
    if missing:
        print(f"WARNING ({tag}): clusters not in ROW_ORDER, appended at end: {sorted(missing)}")
        order += sorted(missing)
    Mo = M.loc[order]

    sns.heatmap(Mo, ax=ax, cmap="Oranges", vmin=0, vmax=vmax, square=True,
               linewidths=0.4, linecolor="white", cbar_kws={"label": "reciprocity-weighted score"},
               yticklabels=list(Mo.index))
    for tick, c in zip(ax.get_yticklabels(), Mo.index):
        tick.set_color(POSITION_COLORS[CLUSTER_POSITION[c]])
    ax.set_xlabel("k.anchor  (weak -> strong integration)", fontsize=LABEL_PT)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=8)
    ax.tick_params(axis="both", labelsize=LABEL_PT, length=2, pad=1.5)
    ax.figure.axes[-1].tick_params(labelsize=LABEL_PT)   # this panel's own colorbar axis

present = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in ROW_ORDER}]
handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=len(handles),
          frameon=False, fontsize=8, title="row label colour  (anatomical position)", title_fontsize=8)
obs_min = min(M.values.min() for M, _ in loaded.values())
fig.suptitle("Per-cluster response to k.anchor, reciprocity-weighted RPCA (hybrid labels)\n"
            "Glut(finch) x Excitatory(chicken)  --  rows in fixed order: song nuclei "
            "(HVCra, HVCx, LMANco, LMANsh, RA), then non-song DACH2, then non-song CACNA1H\n"
            f"colour scale fixed 0-{vmax:.2f} across both panels for comparability "
            f"(observed range {obs_min:.2f}-{vmax:.2f})",
            fontsize=8, y=1.22)
fig.tight_layout()
fig.subplots_adjust(wspace=0.9)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"sweep_heatmap_reciprocal.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT}/sweep_heatmap_reciprocal.pdf / .png  ({fig_w:.1f}x{fig_h:.1f} in)")
