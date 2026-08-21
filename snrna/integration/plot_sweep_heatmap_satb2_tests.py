"""Cluster x k.anchor heatmaps for the two SATB2/mesopallial follow-up analyses, testing
whether Glut-DACH2-HVCra's SCT reciprocal match to the SATB2-family chicken clusters is
real or an artifact of Glut-SATB2-1's exclusion from the main glut-only hybrid subset
(Glut-SATB2-1 is 100% cell-identical to the OLD scheme's Glut-NR-5, confirmed by barcode
crosswalk, which strongly and specifically matched the same chicken SATB2-family types
when present in that analysis).

  1. gg_glutonly_hybrid_plusSATB2: standard glut-only hybrid finch subset + Glut-SATB2-1
     added back in (chicken side unchanged) -- an extra row here.
  2. gg_glutonly_hybrid_noMeso: standard finch subset, chicken reference with
     Ex_SATB2*/Ex_KIAA1217* mesopallial types removed entirely -- same 17-cluster row
     set as the main analysis (only the chicken column pool changed, which doesn't
     affect this per-cluster-best-score view).

Same conventions as plot_sweep_heatmap_glutonly_hybrid_sct_kf200.py (square cells, 6pt
labels, anatomical-position row colours, shared vmax per figure); k.filter is cosmetic
for SCT (see R/cca_transfer.R) so only one k.filter value is plotted per tag.
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
R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/gg_glutonly_hybrid")

BASE_ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]
# Glut-SATB2-1 is not in position_palette.py (excluded from the main analysis) --
# 100% identical cells to the old scheme's Glut-NR-5 (barcode-crosswalk confirmed),
# which was classified non-song ("NR") there, so treated as non-song / "nr"-adjacent
# here too. Given its own distinct colour so it isn't mistaken for an NC/NR row.
EXTRA_POSITION = {"Glut-SATB2-1": "satb2"}
EXTRA_COLOR = {"satb2": "#6a3d9a"}   # common_aesthetics.R's 'ov' slot -- distinct purple, unused elsewhere here
EXTRA_LABEL = {"satb2": "SATB2-1 (=old NR-5)"}


def cluster_pos(c):
    if c in CLUSTER_POSITION:
        return CLUSTER_POSITION[c]
    return EXTRA_POSITION[c]


def pos_color(p):
    return POSITION_COLORS.get(p, EXTRA_COLOR.get(p))


def pos_label(p):
    return POSITION_LABELS.get(p, EXTRA_LABEL.get(p))


def make_figure(tag, row_order, out_stub, title_suffix):
    rows = {}
    for ka in KA:
        f = R / f"{tag}_rpca_SCT_finch_from_mouse_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
        rv = R / f"{tag}_rpca_SCT_mouse_from_finch_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
        if not (f.exists() and rv.exists()):
            print(f"WARNING ({tag}): missing ka={ka} files, skipping")
            continue
        best_score, _, _ = best_per_cluster(f, rv)
        rows[ka] = best_score
    M = pd.DataFrame(rows)
    M = M.drop(index=[c for c in EXCLUDED_CLUSTERS if c in M.index])
    present_ka = list(M.columns)

    order = [c for c in row_order if c in M.index]
    missing = set(M.index) - set(order)
    if missing:
        print(f"WARNING ({tag}): clusters not in row_order, appended at end: {sorted(missing)}")
        order += sorted(missing)
    Mo = M.loc[order, present_ka]

    LABEL_PT = 6
    n_rows = len(order)
    row_pitch_in = LABEL_PT / 72.0 * 1.15
    core_h = row_pitch_in * n_rows
    n_cols = len(present_ka)
    core_w = row_pitch_in * n_cols
    fig_w = core_w + 2.2
    fig_h = core_h + 1.55

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(Mo, ax=ax, cmap="Oranges", vmin=0, vmax=Mo.values.max(), square=True,
               linewidths=0.4, linecolor="white", cbar_kws={"label": "reciprocity-weighted score"},
               yticklabels=list(Mo.index))
    for tick, c in zip(ax.get_yticklabels(), Mo.index):
        tick.set_color(pos_color(cluster_pos(c)))
    ax.set_xlabel("k.anchor  (weak -> strong integration)", fontsize=LABEL_PT)
    ax.set_ylabel("")
    ax.set_title(f"RPCA (SCTransform) -- {title_suffix}", fontsize=8)
    ax.tick_params(axis="both", labelsize=LABEL_PT, length=2, pad=1.5)
    ax.figure.axes[-1].tick_params(labelsize=LABEL_PT)

    present = [p for p in list(POSITION_ORDER) + ["satb2"] if p in {cluster_pos(c) for c in order}]
    handles = [Patch(facecolor=pos_color(p), edgecolor="none", label=pos_label(p)) for p in present]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=min(len(handles), 4),
              frameon=False, fontsize=7.5, title="row label colour", title_fontsize=8)
    fig.suptitle(f"Per-cluster response to k.anchor, reciprocity-weighted RPCA (hybrid, SCT)\n{title_suffix}",
                fontsize=9, y=1.22)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{out_stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{out_stub}.pdf / .png")

    Mout = M.copy()
    Mout["group"] = [song_group(c) if c in CLUSTER_POSITION else "non-song" for c in Mout.index]
    Mout.to_csv(OUT / f"{out_stub}.csv")
    print(f"wrote {OUT}/{out_stub}.csv")
    return M


PLUS_ROW_ORDER = BASE_ROW_ORDER[:5] + ["Glut-SATB2-1"] + BASE_ROW_ORDER[5:]
make_figure("gg_glutonly_hybrid_plusSATB2", PLUS_ROW_ORDER,
           "sweep_heatmap_SCT_plusSATB2_reciprocal",
           "Glut-SATB2-1 (=old Glut-NR-5) added back to the finch subset")
make_figure("gg_glutonly_hybrid_noMeso", BASE_ROW_ORDER,
           "sweep_heatmap_SCT_noMeso_reciprocal",
           "chicken Ex_SATB2*/Ex_KIAA1217* mesopallial types removed")
