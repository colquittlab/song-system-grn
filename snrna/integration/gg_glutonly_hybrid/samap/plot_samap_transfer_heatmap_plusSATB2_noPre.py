"""Chicken x finch label-transfer heatmap for the SAMap run on the plusSATB2/noPre subset
(datasets/samap_bf-adult_gg-adult-glutonly/results_plusSATB2_noPre/), parallel to
plot_transfer_heatmap_plusSATB2.py's RPCA version -- but with an important difference in
what the underlying score actually is.

samap.analysis.get_mapping_scores' MappingTable is NOT a one-directional vote statistic
the way Seurat's TransferData is (which is why reciprocal_score.py's sqrt(fwd*rev)
symmetrization exists at all for the RPCA side) -- it is already perfectly symmetric by
construction: M[bf_i, gg_j] == M[gg_j, bf_i] to floating-point precision, verified
directly on this run's samap_mapping_table.csv (max abs diff between the two directional
submatrices = 0.0). SAMap's alignment score is a mutual-nearest-neighbour edge fraction
between the two clusters, which is inherently an undirected/reciprocal quantity -- there
is no separate "reverse direction" to geometric-mean against, so no correction is applied
or needed here; the raw MappingTable IS the reciprocal score.

Same visual conventions as plot_transfer_heatmap_plusSATB2.py (fixed biological row
order, Glut-SATB2-1 in lilac, columns hierarchically clustered, square cells, 6pt labels,
matching fixed-inch geometry) but a THIRD colour (Greens) distinct from both the RPCA
sweep heatmap (Oranges) and the RPCA transfer heatmap (Blues), since this is a different
method's score entirely.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              EXCLUDED_CLUSTERS)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

LILAC = "#cab2d6"
SATB2_CLUSTER = "Glut-SATB2-1"

IN = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/results/plusSATB2_noPre")
OUT = IN

BASE_ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]
ROW_ORDER = BASE_ROW_ORDER[:13] + [SATB2_CLUSTER] + BASE_ROW_ORDER[13:]

# Fixed-inch geometry -- same constants as plot_transfer_heatmap_plusSATB2.py /
# plot_sweep_heatmap_satb2_tests.py / plot_gsi_sct_heatmaps.py /
# plot_composite_confidence_rank_heatmaps.py, standardized across all of this analysis's
# heatmap types (core geometry was already identical; RIGHT_MARGIN_IN/BAR_* now match the
# other bar-panel heatmaps too).
LABEL_PT = 6
ROW_PITCH_IN = LABEL_PT / 72.0 * 1.15
TOP_MARGIN_IN = 0.85
BOTTOM_MARGIN_IN = 1.7
LEFT_MARGIN_IN = 1.45
RIGHT_MARGIN_IN = 2.7
BAR_GAP_IN = 0.15
BAR_W_IN = 0.45


def row_color(c):
    return LILAC if c == SATB2_CLUSTER else POSITION_COLORS[CLUSTER_POSITION[c]]


def main():
    MT = pd.read_csv(IN / "samap_mapping_table.csv", index_col=0)
    bf_rows = [i for i in MT.index if str(i).startswith("bf_")]
    gg_cols = [c for c in MT.columns if str(c).startswith("gg_")]
    A = MT.loc[bf_rows, gg_cols].copy()
    A.index = [i[3:] for i in A.index]
    A.columns = [c[3:] for c in A.columns]

    # Confirm reciprocity (see module docstring) before using the raw matrix as-is.
    B = MT.loc[gg_cols, bf_rows].copy()
    B.index = [i[3:] for i in B.index]; B.columns = [c[3:] for c in B.columns]
    max_asym = (A.values - B.T.loc[A.index, A.columns].values).__abs__().max()
    print(f"max |A - B.T)| (directional asymmetry check): {max_asym:.2e}")

    S = A.drop(index=[c for c in EXCLUDED_CLUSTERS if c in A.index])
    order = [c for c in ROW_ORDER if c in S.index]
    missing = set(S.index) - set(order)
    if missing:
        print(f"WARNING: rows not in ROW_ORDER, appended at end: {sorted(missing)}")
        order += sorted(missing)
    So = S.loc[order]

    col_link = linkage(So.T.values, method="average", metric="euclidean")
    col_order = So.columns[leaves_list(col_link)]
    Mo = So[col_order]

    n_rows, n_cols = Mo.shape
    core_h = ROW_PITCH_IN * n_rows
    core_w = ROW_PITCH_IN * n_cols
    fig_w = LEFT_MARGIN_IN + core_w + BAR_GAP_IN + BAR_W_IN + RIGHT_MARGIN_IN
    fig_h = TOP_MARGIN_IN + core_h + BOTTOM_MARGIN_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_MARGIN_IN / fig_w, BOTTOM_MARGIN_IN / fig_h, core_w / fig_w, core_h / fig_h])
    sns.heatmap(Mo, ax=ax, cmap="Greens", vmin=0, vmax=Mo.values.max(), square=True,
               linewidths=0.3, linecolor="white", cbar=False, xticklabels=True, yticklabels=True)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for tick, c in zip(ax.get_yticklabels(), Mo.index):
        tick.set_fontsize(LABEL_PT)
        tick.set_color(row_color(c))
    for tick in ax.get_xticklabels():
        tick.set_fontsize(LABEL_PT)
        tick.set_rotation(90)
    ax.tick_params(axis="both", length=2, pad=1.5)

    # Max-per-row bar panel, row-aligned with the heatmap via an identical inverted
    # y-axis -- same convention as plot_gsi_sct_heatmaps.py / plot_composite_confidence_rank_heatmaps.py.
    bar_left = LEFT_MARGIN_IN + core_w + BAR_GAP_IN
    ax_bar = fig.add_axes([bar_left / fig_w, BOTTOM_MARGIN_IN / fig_h, BAR_W_IN / fig_w, core_h / fig_h])
    row_max = Mo.max(axis=1)
    y_pos = np.arange(n_rows) + 0.5
    ax_bar.barh(y_pos, row_max.values, height=0.8, color=[row_color(c) for c in Mo.index])
    ax_bar.set_ylim(n_rows, 0)
    ax_bar.set_xlim(0, row_max.max() * 1.15)
    ax_bar.set_yticks([])
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    ax_bar.set_xlabel("max score", fontsize=LABEL_PT)
    ax_bar.tick_params(axis="x", labelsize=LABEL_PT, length=2, pad=1.5)
    ax_bar.tick_params(axis="y", length=0)

    cbar = fig.colorbar(ax.collections[0], ax=ax_bar, shrink=0.45, aspect=12, pad=0.05)
    cbar.set_label("SAMap alignment score\n(reciprocal by construction)", fontsize=7)
    cbar.ax.tick_params(labelsize=LABEL_PT)

    present_pos = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order if c != SATB2_CLUSTER}]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present_pos]
    handles.append(Patch(facecolor=LILAC, edgecolor="none", label="SATB2-1"))
    bar_right_frac = (bar_left + BAR_W_IN) / fig_w
    heat_mid_frac = (BOTTOM_MARGIN_IN + core_h / 2) / fig_h
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(bar_right_frac + 0.22, heat_mid_frac),
              ncol=1, frameon=False, fontsize=6.5, title="row label colour", title_fontsize=7)

    fig.suptitle("Chicken x finch label-transfer matrix, SAMap alignment score (reciprocal by construction)\n"
                "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed",
                fontsize=8)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"samap_transfer_heatmap.{ext}", dpi=220)
    plt.close(fig)
    print(f"wrote {OUT}/samap_transfer_heatmap.pdf / .png")

    best = pd.DataFrame({
        "samap_top": So.idxmax(axis=1),
        "samap_score": So.max(axis=1).round(4),
    })
    best.to_csv(OUT / "samap_reciprocal_best.csv")
    print(f"wrote {OUT}/samap_reciprocal_best.csv")
    pd.set_option("display.width", 120)
    print(best.sort_values("samap_score", ascending=False).to_string())


if __name__ == "__main__":
    main()
