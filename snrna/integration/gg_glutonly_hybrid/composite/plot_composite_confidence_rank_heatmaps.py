"""Confidence and rank-score heatmaps for the plusSATB2_noPre 3-method composite (GSI +
SAMap + CCA/RPCA), in this analysis's own custom heatmap style -- same conventions as
plot_gsi_sct_heatmaps.py (fixed-inch geometry via fig.add_axes, square cells, fixed
biological row order with Glut-DACH2-HVCra-Pre excluded, columns hierarchically clustered
via scipy linkage for a legible order, black border, max-per-row bar panel coloured by
row label, outside-right legend) -- rather than finch-integration-toolkit's
plot_rank_heatmap.py (a general clustermap-based tool this project's other datasets use,
whose internal gridspec layout doesn't accommodate an extra bar panel without disturbing
its shared, carefully-tuned sizing math). No agreement dots here, per explicit request
(this analysis's own confidence/rank-score reads are being examined directly via the
per-method comparison heatmaps already built, not via the dot overlay).

confidence_matrix and rank_score_matrix are both [0,1] magnitudes (not signed), so
sequential Oranges is used (matching plot_rank_heatmap.py's own convention: DIVERGING
only for signed data), with a fixed 0-1 scale for cross-panel comparability between the
two heatmaps.
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

D = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/composite/results/composite_plusSATB2_noPre")
OUT = D

TESTS = [
    ("composite_confidence_matrix.csv", "confidence", "composite confidence (magnitude)"),
    ("composite_rank_score_matrix.csv", "rank_score", "composite rank score (agreement)"),
]

BASE_ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]
ROW_ORDER = BASE_ROW_ORDER[:13] + [SATB2_CLUSTER] + BASE_ROW_ORDER[13:]

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
    for fname, stub, title_suffix in TESTS:
        C = pd.read_csv(D / fname, index_col=0)
        C = C.drop(index=[c for c in EXCLUDED_CLUSTERS if c in C.index])

        order = [c for c in ROW_ORDER if c in C.index]
        missing = set(C.index) - set(order)
        if missing:
            print(f"WARNING {stub}: rows not in ROW_ORDER, appended at end: {sorted(missing)}")
            order += sorted(missing)
        So = C.loc[order]

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
        sns.heatmap(Mo, ax=ax, cmap="Oranges", vmin=0, vmax=1, square=True,
                   linewidths=0.3, linecolor="white", cbar=False,
                   xticklabels=True, yticklabels=True)
        for side in ax.spines.values():
            side.set_visible(True)
            side.set_edgecolor("black")
            side.set_linewidth(1.0)
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
        # y-axis -- same convention as plot_gsi_sct_heatmaps.py.
        bar_left = LEFT_MARGIN_IN + core_w + BAR_GAP_IN
        ax_bar = fig.add_axes([bar_left / fig_w, BOTTOM_MARGIN_IN / fig_h, BAR_W_IN / fig_w, core_h / fig_h])
        row_max = Mo.max(axis=1)
        y_pos = np.arange(n_rows) + 0.5
        ax_bar.barh(y_pos, row_max.values, height=0.8, color=[row_color(c) for c in Mo.index])
        ax_bar.set_ylim(n_rows, 0)
        ax_bar.set_xlim(0, max(1.0, row_max.max() * 1.15))
        ax_bar.set_yticks([])
        for side in ("top", "right", "left"):
            ax_bar.spines[side].set_visible(False)
        ax_bar.set_xlabel(f"max {stub}", fontsize=LABEL_PT)
        ax_bar.tick_params(axis="x", labelsize=LABEL_PT, length=2, pad=1.5)
        ax_bar.tick_params(axis="y", length=0)

        cbar = fig.colorbar(ax.collections[0], ax=ax_bar, shrink=0.45, aspect=12, pad=0.05)
        cbar.set_label(stub.replace("_", " "), fontsize=7)
        cbar.ax.tick_params(labelsize=LABEL_PT)

        present_pos = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order if c != SATB2_CLUSTER}]
        handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present_pos]
        if SATB2_CLUSTER in order:
            handles.append(Patch(facecolor=LILAC, edgecolor="none", label="SATB2-1"))
        bar_right_frac = (bar_left + BAR_W_IN) / fig_w
        heat_mid_frac = (BOTTOM_MARGIN_IN + core_h / 2) / fig_h
        fig.legend(handles=handles, loc="center left", bbox_to_anchor=(bar_right_frac + 0.22, heat_mid_frac),
                  ncol=1, frameon=False, fontsize=6.5, title="row label colour", title_fontsize=7)

        fig.suptitle("plusSATB2_noPre 3-method composite (GSI + SAMap + CCA/RPCA)\n"
                    f"{title_suffix}",
                    fontsize=8)
        for ext in ("pdf", "png"):
            fig.savefig(OUT / f"{stub}_heatmap.{ext}", dpi=220)
        plt.close(fig)
        print(f"wrote {OUT}/{stub}_heatmap.pdf / .png")


if __name__ == "__main__":
    main()
