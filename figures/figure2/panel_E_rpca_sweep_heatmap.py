#!/usr/bin/env python3
"""Figure 2, panel E -- strength of the best finch-chicken match of each finch glutamatergic type across integration strengths.

Seurat RPCA label transfer (SCTransform data, 40 dimensions) between the finch glutamatergic types (Glut-SATB2-1 kept,
Glut-DACH2-HVCra-Int left out) and the chicken excitatory types of Zaremba et al. 2025 (precursor types Ex_Pre_KCNH7 and
Ex_Pre_SATB2 removed), run at k.anchor = 2, 5, 10, 20, 30, 50, 75, 100. For each run and each pair of types, the score in each
direction is the fraction of one type's cells assigned to the other type times their mean prediction score; the reciprocal
score is the geometric mean of the two directions (common.symmetrize). Each cell of the heatmap is a finch type's highest
reciprocal score over all chicken types ("max reciprocal transfer score"). The figure marks k.anchor = 20, the run shown in F.
Same computation and drawing as the plusSATB2_noPre call of snrna/integration/gg_glutonly_hybrid/rpca/plot_sweep_heatmap_satb2_tests.py,
which wrote
  snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid/sweep_heatmap_SCT_plusSATB2_noPre_reciprocal.{pdf,png,csv}
Labels, legends and the k.anchor = 20 box were restyled by hand in the figure.

Input: the forward and reverse label-transfer matrices of each k.anchor (common.matrix_paths).
Output: FIG2_OUT/panel_E_rpca_sweep_heatmap.{pdf,png} and panel_E_rpca_sweep_heatmap.csv (finch type x k.anchor, with song group)

  python figures/figure2/panel_E_rpca_sweep_heatmap.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

# Fixed-inch geometry of the original (rows land at the same height as in panel F's script)
LABEL_PT = 6
ROW_PITCH_IN = LABEL_PT / 72.0 * 1.15
TOP_MARGIN_IN, BOTTOM_MARGIN_IN, LEFT_MARGIN_IN, RIGHT_MARGIN_IN = 0.85, 1.7, 1.45, 2.3


def pos_of(cl):
    return "satb2" if cl == c.SATB2_CLUSTER else c.CLUSTER_POSITION[cl]


def pos_color(p):
    return c.LILAC if p == "satb2" else c.POSITION_COLORS[p]


def pos_label(p):
    return "SATB2-1 (=old NR-5)" if p == "satb2" else c.POSITION_LABELS[p]


rows = {}
for ka in c.KA_VALUES:
    f, r = c.matrix_paths(ka)
    best, _, _ = c.best_per_cluster(f, r)
    rows[ka] = best
M = pd.DataFrame(rows)
M = M.drop(index=[cl for cl in c.EXCLUDED_CLUSTERS if cl in M.index])
order = [cl for cl in c.ROW_ORDER if cl in M.index]
assert set(order) == set(M.index), f"finch types not in ROW_ORDER: {sorted(set(M.index) - set(order))}"
Mo = M.loc[order, c.KA_VALUES]

n_rows, n_cols = Mo.shape
core_h, core_w = ROW_PITCH_IN * n_rows, ROW_PITCH_IN * n_cols
fig_w = LEFT_MARGIN_IN + core_w + RIGHT_MARGIN_IN
fig_h = TOP_MARGIN_IN + core_h + BOTTOM_MARGIN_IN
fig = plt.figure(figsize=(fig_w, fig_h))
ax = fig.add_axes([LEFT_MARGIN_IN / fig_w, BOTTOM_MARGIN_IN / fig_h, core_w / fig_w, core_h / fig_h])
sns.heatmap(Mo, ax=ax, cmap="Oranges", vmin=0, vmax=Mo.values.max(), square=True,
            linewidths=0.4, linecolor="white", cbar=False, yticklabels=list(Mo.index))
cbar = fig.colorbar(ax.collections[0], ax=ax, shrink=0.45, aspect=12, pad=0.03)
cbar.set_label("reciprocity-weighted score", fontsize=7)
cbar.ax.tick_params(labelsize=LABEL_PT)
for tick, cl in zip(ax.get_yticklabels(), Mo.index):
    tick.set_fontsize(LABEL_PT)
    tick.set_color(pos_color(pos_of(cl)))
ax.set_xlabel("k.anchor  (weak -> strong integration)", fontsize=LABEL_PT)
ax.set_ylabel("")
ax.tick_params(axis="both", labelsize=LABEL_PT, length=2, pad=1.5)

present = [p for p in c.POSITION_ORDER + ["satb2"] if p in {pos_of(cl) for cl in order}]
handles = [Patch(facecolor=pos_color(p), edgecolor="none", label=pos_label(p)) for p in present]
fig.legend(handles=handles, loc="center left",
           bbox_to_anchor=((LEFT_MARGIN_IN + core_w) / fig_w + 0.16, (BOTTOM_MARGIN_IN + core_h / 2) / fig_h),
           ncol=1, frameon=False, fontsize=7, title="row label colour", title_fontsize=7.5)
fig.suptitle("Per-cluster response to k.anchor, reciprocity-weighted RPCA (SCTransform)\n"
             "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed\n"
             f"k.filter={c.KF} (nominal only, forced NA under SCT), dims={c.DIMS}",
             fontsize=8)
out = c.FIG2_OUT / "panel_E_rpca_sweep_heatmap"
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}", dpi=220)
plt.close(fig)

Mout = M.copy()
Mout["group"] = [c.song_group(cl) for cl in Mout.index]
Mout.to_csv(f"{out}.csv")
print(f"wrote {out}.pdf / .png / .csv")
