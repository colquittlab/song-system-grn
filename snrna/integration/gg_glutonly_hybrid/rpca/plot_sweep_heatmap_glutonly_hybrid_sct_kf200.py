"""Cluster x k.anchor heatmap of the reciprocity-weighted RPCA sweep, SCTransform variant
at k.filter=200 -- companion to plot_sweep_heatmap_glutonly_hybrid_sct.py's k.filter=50
("manuscript config") version, added per explicit request to make the SCT track directly
comparable to the log-norm analysis (which used k.filter=200 throughout). Both k.filter
values are kept side by side under SCT, not one replacing the other. Same conventions
otherwise (fixed biological row order, Glut-DACH2-HVCra-Pre excluded, square cells, 6pt
labels, anatomical-position row colours).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

import sys
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
from reciprocal_score import best_per_cluster
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              EXCLUDED_CLUSTERS, song_group)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
KA = [2, 5, 10, 20, 30, 50, 75, 100]
KF = 200
DIMS = 40
R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid")
TAG = "gg_glutonly_hybrid"

ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]

rows = {}
for ka in KA:
    f = R / f"{TAG}_rpca_SCT_finch_from_mouse_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
    rv = R / f"{TAG}_rpca_SCT_mouse_from_finch_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
    if not (f.exists() and rv.exists()):
        print(f"WARNING: missing ka={ka} files, skipping")
        continue
    best_score, _, _ = best_per_cluster(f, rv)
    rows[ka] = best_score
M = pd.DataFrame(rows)
M = M.drop(index=[c for c in EXCLUDED_CLUSTERS if c in M.index])
present_ka = list(M.columns)

order = [c for c in ROW_ORDER if c in M.index]
missing = set(M.index) - set(order)
if missing:
    print(f"WARNING: clusters not in ROW_ORDER, appended at end: {sorted(missing)}")
    order += sorted(missing)
Mo = M.loc[order]

LABEL_PT = 6
n_rows = len(order)
row_pitch_in = LABEL_PT / 72.0 * 1.15
core_h = row_pitch_in * n_rows
n_cols = len(present_ka)
core_w = row_pitch_in * n_cols
fig_w = core_w + 1.9
fig_h = core_h + 1.55

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
sns.heatmap(Mo, ax=ax, cmap="Oranges", vmin=0, vmax=Mo.values.max(), square=True,
           linewidths=0.4, linecolor="white", cbar_kws={"label": "reciprocity-weighted score"},
           yticklabels=list(Mo.index))
for tick, c in zip(ax.get_yticklabels(), Mo.index):
    tick.set_color(POSITION_COLORS[CLUSTER_POSITION[c]])
ax.set_xlabel("k.anchor  (weak -> strong integration)", fontsize=LABEL_PT)
ax.set_ylabel("")
ax.set_title(f"RPCA (SCTransform, k.filter={KF} -- matches log-norm's k.filter, dims={DIMS})", fontsize=8)
ax.tick_params(axis="both", labelsize=LABEL_PT, length=2, pad=1.5)
ax.figure.axes[-1].tick_params(labelsize=LABEL_PT)

present = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order}]
handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=len(handles),
          frameon=False, fontsize=8, title="row label colour  (anatomical position)", title_fontsize=8)
fig.suptitle("Per-cluster response to k.anchor, reciprocity-weighted RPCA -- SCTransform, k.filter=200\n"
            "Glut(finch) x Excitatory(chicken)", fontsize=9, y=1.20)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"sweep_heatmap_SCT_kf200_reciprocal.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT}/sweep_heatmap_SCT_kf200_reciprocal.pdf / .png")

stats_rows = []
for ka in present_ka:
    song_vals = M.loc[[c for c in M.index if song_group(c) == "song"], ka]
    nonsong_vals = M.loc[[c for c in M.index if song_group(c) == "non-song"], ka]
    u, p = mannwhitneyu(song_vals, nonsong_vals, alternative="less")
    stats_rows.append({"k_anchor": ka, "song_mean": song_vals.mean(), "nonsong_mean": nonsong_vals.mean(),
                       "U": u, "p_one_sided": p})
stats = pd.DataFrame(stats_rows)
stats.to_csv(OUT / "sweep_percluster_SCT_kf200_reciprocal_stats.csv", index=False)
Mout = M.copy()
Mout["group"] = [song_group(c) for c in Mout.index]
Mout.to_csv(OUT / "sweep_percluster_SCT_kf200_reciprocal.csv")
pd.set_option("display.width", 120)
print(stats.to_string(index=False))
