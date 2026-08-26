"""Full chicken x finch label-transfer heatmaps (reciprocity-weighted RPCA scores) for
every k.anchor in the plusSATB2 SCT sweep and its Ex_Pre_*-excluded follow-up -- unlike
plot_sweep_heatmap_satb2_tests.py (which collapses each finch row to its single
best-scoring chicken column, to compare across k.anchor), this shows the COMPLETE
finch-row x chicken-column matrix at each individual k.anchor, so partial/competing
matches are visible too.

  1. gg_glutonly_hybrid_plusSATB2: standard glut-only hybrid finch subset + Glut-SATB2-1
     added back in; chicken side unchanged (gg_adult_ex.h5ad).
  2. gg_glutonly_hybrid_plusSATB2_noPre: same finch subset as (1); chicken reference with
     BOTH "precursor" populations removed (Ex_Pre_KCNH7 -- noticed as an
     unexpectedly-present population in (1)'s own heatmap -- and Ex_Pre_SATB2, caught on a
     second pass after Ex_Pre_KCNH7 alone turned out to have been an incomplete fix: any
     Ex_Pre_* cluster is arguably not a definitive excitatory reference type).

S[i,j] = sqrt(fwd[i,j] * rev[j,i]) (reciprocal_score.symmetrize) -- a pair only scores
high if both label-transfer directions independently nominate it.

Row order is fixed biological order (Glut-SATB2-1 inserted after the DACH2 clusters,
before the CACNA1H clusters); columns are hierarchically clustered (via scipy linkage, no
dendrogram drawn -- only used to pick a legible column order) so similar chicken columns
group together without disturbing the fixed row order. Only Glut-SATB2-1's row label is
coloured light purple/lilac (#cab2d6, common_aesthetics.R's 'fieldl' slot) -- chicken
column labels are left uncoloured. Square cells, 6pt labels, compressed spacing.

Geometry (LABEL_PT/ROW_PITCH_IN and the four *_MARGIN_IN constants) is fixed-inch and
duplicated verbatim in plot_sweep_heatmap_satb2_tests.py, deliberately, so that for the
same finch row set the heatmap CORE (row 1's top edge, row pitch, row N's bottom edge)
lands at identical absolute figure positions in both scripts' output -- letting the two
be stacked/overlaid with rows aligned. Axes are placed via fig.add_axes at an explicit
inch-derived rect rather than tight_layout/constrained_layout, and savefig does NOT use
bbox_inches="tight" -- either would re-crop the canvas to content and silently break the
alignment this buys.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

import sys
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
from reciprocal_score import symmetrize
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              EXCLUDED_CLUSTERS)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

LILAC = "#cab2d6"
SATB2_CLUSTER = "Glut-SATB2-1"

KA_VALUES = [2, 5, 10, 20, 30, 50, 75, 100]
KF, DIMS = 200, 40
R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid")

BASE_ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]
ROW_ORDER = BASE_ROW_ORDER[:13] + [SATB2_CLUSTER] + BASE_ROW_ORDER[13:]

# Shared fixed-inch geometry -- kept identical to plot_sweep_heatmap_satb2_tests.py so
# the two heatmap types' rows land at the same absolute position (see module docstring).
LABEL_PT = 6
ROW_PITCH_IN = LABEL_PT / 72.0 * 1.15
TOP_MARGIN_IN = 0.85
BOTTOM_MARGIN_IN = 1.7
LEFT_MARGIN_IN = 1.45
RIGHT_MARGIN_IN = 1.7

TESTS = [
    ("gg_glutonly_hybrid_plusSATB2", "plusSATB2",
     "Glut-SATB2-1 (=old Glut-NR-5) added back to the finch subset"),
    ("gg_glutonly_hybrid_plusSATB2_noPre", "plusSATB2_noPre",
     "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed"),
]


def row_color(c):
    return LILAC if c == SATB2_CLUSTER else POSITION_COLORS[CLUSTER_POSITION[c]]


def make_figure(tag, stub_tag, title_suffix, ka):
    f = R / f"{tag}_rpca_SCT_finch_from_mouse_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
    rv = R / f"{tag}_rpca_SCT_mouse_from_finch_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
    if not (f.exists() and rv.exists()):
        print(f"SKIP {tag} ka={ka}: matrix files not found")
        return
    S = symmetrize(f, rv)
    S = S.drop(index=[c for c in EXCLUDED_CLUSTERS if c in S.index])

    order = [c for c in ROW_ORDER if c in S.index]
    missing = set(S.index) - set(order)
    if missing:
        print(f"WARNING {tag} ka={ka}: rows not in ROW_ORDER, appended at end: {sorted(missing)}")
        order += sorted(missing)
    So = S.loc[order]

    col_link = linkage(So.T.values, method="average", metric="euclidean")
    col_order = So.columns[leaves_list(col_link)]
    Mo = So[col_order]

    n_rows, n_cols = Mo.shape
    core_h = ROW_PITCH_IN * n_rows
    core_w = ROW_PITCH_IN * n_cols
    fig_w = LEFT_MARGIN_IN + core_w + RIGHT_MARGIN_IN
    fig_h = TOP_MARGIN_IN + core_h + BOTTOM_MARGIN_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_MARGIN_IN / fig_w, BOTTOM_MARGIN_IN / fig_h, core_w / fig_w, core_h / fig_h])
    sns.heatmap(Mo, ax=ax, cmap="Blues", vmin=0, vmax=Mo.values.max(), square=True,
               linewidths=0.3, linecolor="white", cbar=False,
               xticklabels=True, yticklabels=True)
    cbar = fig.colorbar(ax.collections[0], ax=ax, shrink=0.45, aspect=12, pad=0.03)
    cbar.set_label("reciprocity-weighted score", fontsize=7)
    cbar.ax.tick_params(labelsize=LABEL_PT)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for tick, c in zip(ax.get_yticklabels(), Mo.index):
        tick.set_fontsize(LABEL_PT)
        tick.set_color(row_color(c))
    for tick in ax.get_xticklabels():
        tick.set_fontsize(LABEL_PT)
        tick.set_rotation(90)
    ax.tick_params(axis="both", length=2, pad=1.5)

    present_pos = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order if c != SATB2_CLUSTER}]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present_pos]
    handles.append(Patch(facecolor=LILAC, edgecolor="none", label="SATB2-1"))
    heat_right_frac = (LEFT_MARGIN_IN + core_w) / fig_w
    heat_mid_frac = (BOTTOM_MARGIN_IN + core_h / 2) / fig_h
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(heat_right_frac + 0.16, heat_mid_frac),
              ncol=1, frameon=False, fontsize=6.5, title="row label colour", title_fontsize=7)

    fig.suptitle("Chicken x finch label-transfer matrix, reciprocity-weighted RPCA (SCTransform)\n"
                f"{title_suffix}\n"
                f"k.anchor={ka}, k.filter={KF} (nominal only, forced NA under SCT), dims={DIMS}",
                fontsize=8)
    stub = f"transfer_heatmap_{stub_tag}_ka{ka}_kf{KF}_d{DIMS}_SCT"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220)
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


if __name__ == "__main__":
    for tag, stub_tag, title_suffix in TESTS:
        for ka in KA_VALUES:
            make_figure(tag, stub_tag, title_suffix, ka)
