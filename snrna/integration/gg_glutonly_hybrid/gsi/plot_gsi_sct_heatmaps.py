"""Full chicken x finch GSI-SCT correlation heatmaps, one per finch/chicken subset
combination (base, plusSATB2, noMeso, plusSATB2_noPre) -- same visual conventions as
plot_transfer_heatmap_plusSATB2.py's RPCA transfer-matrix heatmaps (fixed-inch geometry
via fig.add_axes, no dendrogram/tight_layout, square cells, 6pt labels, fixed biological
row order with Glut-SATB2-1 after the DACH2 clusters coloured lilac, columns
hierarchically clustered via scipy linkage for a legible order, outside-right legend,
shrunk colourbar).

Unlike the RPCA/SAMap transfer heatmaps, GSI's matrix is a Spearman CORRELATION (signed,
can be negative) rather than a bounded [0,1] magnitude -- no reciprocity symmetrization
is needed (unlike Seurat's one-directional TransferData, GSI's cross-species correlation
is computed directly from both species' specificity profiles at once) but the colour
ramp must be DIVERGING through white, not sequential from 0, per this project's own
established convention for signed data (see finch-integration-toolkit/plot_rank_heatmap.py:
"a sequential ramp on signed data makes a strong negative and a weak positive look
alike"). Colour scale is shared (symmetric vmin/vmax from the global max |r| across all
four subsets) so the four separate figures stay directly comparable to each other.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              EXCLUDED_CLUSTERS)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

# Established project convention for signed/correlation data (finch-integration-toolkit's
# plot_rank_heatmap.py) -- reused verbatim rather than inventing a new hue.
DIVERGING = LinearSegmentedColormap.from_list("corr_div", ["#2a78d6", "#ffffff", "#eb6834"])

LILAC = "#cab2d6"
SATB2_CLUSTER = "Glut-SATB2-1"

R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/gsi/results")
OUT = R

TESTS = [
    ("base", "gsi_corr_gg_glutonly_hybrid_SCT.csv", "base (standard glut-only hybrid subset)"),
    ("plusSATB2", "gsi_corr_gg_glutonly_hybrid_plusSATB2_SCT.csv", "Glut-SATB2-1 added back to the finch subset"),
    ("noMeso", "gsi_corr_gg_glutonly_hybrid_noMeso_SCT.csv", "chicken Ex_SATB2*/Ex_KIAA1217* mesopallial types removed"),
    ("plusSATB2_noPre", "gsi_corr_gg_glutonly_hybrid_plusSATB2_noPre_SCT.csv",
     "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed"),
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
# Max-per-row bar panel, inserted between the heatmap and the colourbar/legend.
BAR_GAP_IN = 0.15
BAR_W_IN = 0.45


def row_color(c):
    return LILAC if c == SATB2_CLUSTER else POSITION_COLORS[CLUSTER_POSITION[c]]


def main():
    mats = {}
    for tag, fname, _ in TESTS:
        C = pd.read_csv(R / fname, index_col=0)
        C = C.drop(index=[c for c in EXCLUDED_CLUSTERS if c in C.index])
        mats[tag] = C
    vmax = max(C.values.max() for C in mats.values())
    vmin = min(C.values.min() for C in mats.values())
    absmax = max(abs(vmin), abs(vmax))
    print(f"shared color scale: -{absmax:.3f} to +{absmax:.3f}")

    for tag, fname, title_suffix in TESTS:
        C = mats[tag]
        order = [c for c in ROW_ORDER if c in C.index]
        missing = set(C.index) - set(order)
        if missing:
            print(f"WARNING {tag}: rows not in ROW_ORDER, appended at end: {sorted(missing)}")
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
        sns.heatmap(Mo, ax=ax, cmap=DIVERGING, vmin=-absmax, vmax=absmax, square=True,
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

        # Max-per-row bar panel: each finch cluster's single best (max) GSI r across all
        # chicken clusters in this subset, coloured by the same row label colour as the
        # heatmap -- aligned row-for-row via an identical (n_rows, 0) inverted y-axis.
        bar_left = LEFT_MARGIN_IN + core_w + BAR_GAP_IN
        ax_bar = fig.add_axes([bar_left / fig_w, BOTTOM_MARGIN_IN / fig_h, BAR_W_IN / fig_w, core_h / fig_h])
        row_max = Mo.max(axis=1)
        y_pos = np.arange(n_rows) + 0.5
        ax_bar.barh(y_pos, row_max.values, height=0.8, color=[row_color(c) for c in Mo.index])
        ax_bar.set_ylim(n_rows, 0)
        bar_lo = min(0.0, row_max.min() * 1.1)
        bar_hi = row_max.max() * 1.15
        ax_bar.set_xlim(bar_lo, bar_hi)
        if bar_lo < 0:
            ax_bar.axvline(0, color="black", linewidth=0.5)
        ax_bar.set_yticks([])
        for side in ("top", "right", "left"):
            ax_bar.spines[side].set_visible(False)
        ax_bar.set_xlabel("max r", fontsize=LABEL_PT)
        ax_bar.tick_params(axis="x", labelsize=LABEL_PT, length=2, pad=1.5)
        ax_bar.tick_params(axis="y", length=0)

        cbar = fig.colorbar(ax.collections[0], ax=ax_bar, shrink=0.45, aspect=12, pad=0.05)
        cbar.set_label("GSI Spearman r", fontsize=7)
        cbar.ax.tick_params(labelsize=LABEL_PT)

        present_pos = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order if c != SATB2_CLUSTER}]
        handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present_pos]
        if SATB2_CLUSTER in order:
            handles.append(Patch(facecolor=LILAC, edgecolor="none", label="SATB2-1"))
        bar_right_frac = (bar_left + BAR_W_IN) / fig_w
        heat_mid_frac = (BOTTOM_MARGIN_IN + core_h / 2) / fig_h
        fig.legend(handles=handles, loc="center left", bbox_to_anchor=(bar_right_frac + 0.22, heat_mid_frac),
                  ncol=1, frameon=False, fontsize=6.5, title="row label colour", title_fontsize=7)

        fig.suptitle("Chicken x finch GSI-SCT correlation matrix\n"
                    f"{title_suffix}",
                    fontsize=8)
        stub = f"gsi_sct_heatmap_{tag}"
        for ext in ("pdf", "png"):
            fig.savefig(OUT / f"{stub}.{ext}", dpi=220)
        plt.close(fig)
        print(f"wrote {OUT}/{stub}.pdf / .png")


if __name__ == "__main__":
    main()
