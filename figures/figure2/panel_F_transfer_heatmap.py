#!/usr/bin/env python3
"""Figure 2, panel F -- reciprocal transfer score between every finch glutamatergic type and every chicken excitatory type at k.anchor = 20.

The full finch-type x chicken-type matrix of reciprocal scores (geometric mean of the two label-transfer directions,
common.symmetrize) from the RPCA run of panel E at k.anchor = 20 (SCTransform data, 40 dimensions; finch Glut-SATB2-1 kept,
Glut-DACH2-HVCra-Int left out; chicken Ex_Pre_KCNH7 and Ex_Pre_SATB2 removed). Rows in panel E's order; chicken columns
grouped by the marker in their name (DACH2, TSHZ2, CACNA1H, SATB2, KIAA1217, BCL6; alphabetical within a group,
common.chicken_col_order). Colour scale from 0 to the matrix maximum.
Same computation and drawing as the plusSATB2_noPre, k.anchor = 20 call of
snrna/integration/gg_glutonly_hybrid/rpca/plot_transfer_heatmap_plusSATB2.py -- the version with the fixed chicken-column order
(main checkout working copy, Sep 25 2026; the committed version orders the columns by hierarchical clustering instead), which wrote
  snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid/transfer_heatmap_plusSATB2_noPre_ka20_kf200_d40_SCT.{pdf,png}
In the figure the finch row labels are dropped (they are shared with panel E) and the legend restyled.

Input: the forward and reverse label-transfer matrices at k.anchor = 20 (common.matrix_paths).
Output: FIG2_OUT/panel_F_transfer_heatmap_ka20.{pdf,png} and panel_F_transfer_heatmap_ka20.csv (the plotted matrix)

  python figures/figure2/panel_F_transfer_heatmap.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
KA = 20

# Fixed-inch geometry of the original (rows land at the same height as in panel E's script)
LABEL_PT = 6
ROW_PITCH_IN = LABEL_PT / 72.0 * 1.15
TOP_MARGIN_IN, BOTTOM_MARGIN_IN, LEFT_MARGIN_IN, RIGHT_MARGIN_IN = 0.85, 1.7, 1.45, 1.7

f, r = c.matrix_paths(KA)
S = c.symmetrize(f, r)
S = S.drop(index=[cl for cl in c.EXCLUDED_CLUSTERS if cl in S.index])
order = [cl for cl in c.ROW_ORDER if cl in S.index]
assert set(order) == set(S.index), f"finch types not in ROW_ORDER: {sorted(set(S.index) - set(order))}"
Mo = S.loc[order, c.chicken_col_order(S.columns)]

n_rows, n_cols = Mo.shape
core_h, core_w = ROW_PITCH_IN * n_rows, ROW_PITCH_IN * n_cols
fig_w = LEFT_MARGIN_IN + core_w + RIGHT_MARGIN_IN
fig_h = TOP_MARGIN_IN + core_h + BOTTOM_MARGIN_IN
fig = plt.figure(figsize=(fig_w, fig_h))
ax = fig.add_axes([LEFT_MARGIN_IN / fig_w, BOTTOM_MARGIN_IN / fig_h, core_w / fig_w, core_h / fig_h])
sns.heatmap(Mo, ax=ax, cmap="Blues", vmin=0, vmax=Mo.values.max(), square=True,
            linewidths=0.3, linecolor="white", cbar=False, xticklabels=True, yticklabels=True)
cbar = fig.colorbar(ax.collections[0], ax=ax, shrink=0.45, aspect=12, pad=0.03)
cbar.set_label("reciprocity-weighted score", fontsize=7)
cbar.ax.tick_params(labelsize=LABEL_PT)
ax.set_xlabel("")
ax.set_ylabel("")
for tick, cl in zip(ax.get_yticklabels(), Mo.index):
    tick.set_fontsize(LABEL_PT)
    tick.set_color(c.row_color(cl))
for tick in ax.get_xticklabels():
    tick.set_fontsize(LABEL_PT)
    tick.set_rotation(90)
ax.tick_params(axis="both", length=2, pad=1.5)

present_pos = [p for p in c.POSITION_ORDER if p in {c.CLUSTER_POSITION[cl] for cl in order if cl != c.SATB2_CLUSTER}]
handles = [Patch(facecolor=c.POSITION_COLORS[p], edgecolor="none", label=c.POSITION_LABELS[p]) for p in present_pos]
handles.append(Patch(facecolor=c.LILAC, edgecolor="none", label="SATB2-1"))
fig.legend(handles=handles, loc="center left",
           bbox_to_anchor=((LEFT_MARGIN_IN + core_w) / fig_w + 0.16, (BOTTOM_MARGIN_IN + core_h / 2) / fig_h),
           ncol=1, frameon=False, fontsize=6.5, title="row label colour", title_fontsize=7)
fig.suptitle("Chicken x finch label-transfer matrix, reciprocity-weighted RPCA (SCTransform)\n"
             "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed\n"
             f"k.anchor={KA}, k.filter={c.KF} (nominal only, forced NA under SCT), dims={c.DIMS}",
             fontsize=8)
out = c.FIG2_OUT / f"panel_F_transfer_heatmap_ka{KA}"
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}", dpi=220)
plt.close(fig)
Mo.to_csv(f"{out}.csv")
print(f"wrote {out}.pdf / .png / .csv")
