"""Summary heatmap comparing GSI-SCT top matches across the four finch/chicken subset
combinations now used throughout this analysis (base, plusSATB2, noMeso, plusSATB2_noPre)
-- a direct visualization of the composition-dependency question: for each finch cluster,
does its single best-matching chicken cluster (and the strength of that match) change
depending on which cell types are included in the comparison?

Same visual conventions as plot_transfer_heatmap_plusSATB2.py / plot_sweep_heatmap_satb2_tests.py
(fixed-inch geometry via fig.add_axes, no dendrogram/tight_layout, square cells, 6pt
labels, position-palette row colours, "outside right" legend, shrunk colourbar) -- but
columns here are the four SUBSETS (not chicken clusters), and each cell is annotated with
which chicken cluster matched (sns.heatmap's annot=, not a separate axis), since the
whole point of this comparison is identity, not just magnitude.

Cell colour: GSI Spearman correlation (r) of that finch cluster's best chicken match in
that subset -- sequential Purples, a fourth colour distinct from the RPCA sweep heatmap
(Oranges), RPCA transfer heatmap (Blues), and SAMap transfer heatmap (Greens).
Rows: finch clusters (fixed biological order, Glut-SATB2-1 after the DACH2 clusters --
present only in plusSATB2/plusSATB2_noPre, since those are the only two subsets keeping
that finch cluster; base/noMeso get a blank cell for that row).
Cells whose matched identity differs from the finch cluster's own MODE (most common match
across the subsets it appears in) get a bold red border -- a direct visual flag for
composition-sensitive calls, the finding this comparison exists to surface.
"""
from pathlib import Path
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd
import numpy as np
import seaborn as sns

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              EXCLUDED_CLUSTERS)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

LILAC = "#cab2d6"
SATB2_CLUSTER = "Glut-SATB2-1"

R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/gsi/results")
OUT = R

TESTS = [
    ("base", "gsi_corr_gg_glutonly_hybrid_SCT.csv", "base"),
    ("plusSATB2", "gsi_corr_gg_glutonly_hybrid_plusSATB2_SCT.csv", "+SATB2-1"),
    ("noMeso", "gsi_corr_gg_glutonly_hybrid_noMeso_SCT.csv", "noMeso"),
    ("plusSATB2_noPre", "gsi_corr_gg_glutonly_hybrid_plusSATB2_noPre_SCT.csv", "+SATB2-1,\nnoPre"),
]

BASE_ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]
ROW_ORDER = BASE_ROW_ORDER[:13] + [SATB2_CLUSTER] + BASE_ROW_ORDER[13:]

LABEL_PT = 6
# Larger than the other heatmaps' 6pt-based ROW_PITCH_IN (~0.10in) -- those never put
# text INSIDE cells (colour + outside tick labels only); here every cell carries a
# two-line annotation (chicken cluster name + r), which needs real room to stay legible.
ROW_PITCH_IN = 0.72
TOP_MARGIN_IN = 0.85
BOTTOM_MARGIN_IN = 0.55
LEFT_MARGIN_IN = 1.45
RIGHT_MARGIN_IN = 1.9


def row_color(c):
    return LILAC if c == SATB2_CLUSTER else POSITION_COLORS[CLUSTER_POSITION[c]]


def short(c):
    return c.replace("Ex_", "")


def wrap_id(cid):
    """Break a chicken cluster's short name at its first underscore onto two lines
    (e.g. "DACH2_MGAT4C" -> "DACH2\nMGAT4C") so it fits the cell width at readable
    font size -- a single line at 6pt overflows most >8-character names."""
    if "_" in cid and len(cid) > 8:
        fam, gene = cid.split("_", 1)
        return f"{fam}\n{gene}"
    return cid


def main():
    top = {}   # {finch_cluster: {tag: (chicken_short, r)}}
    for tag, fname, _ in TESTS:
        C = pd.read_csv(R / fname, index_col=0)
        for fc in C.index:
            if fc in EXCLUDED_CLUSTERS:
                continue
            row = C.loc[fc]
            best_col = row.idxmax()
            top.setdefault(fc, {})[tag] = (short(best_col), round(float(row[best_col]), 3))

    order = [c for c in ROW_ORDER if c in top]
    missing = set(top) - set(order)
    if missing:
        print(f"WARNING: clusters not in ROW_ORDER, appended at end: {sorted(missing)}")
        order += sorted(missing)

    tags = [t for t, _, _ in TESTS]
    col_labels = [lbl for _, _, lbl in TESTS]
    n_rows, n_cols = len(order), len(tags)

    R_mat = pd.DataFrame(np.nan, index=order, columns=tags)
    annot = pd.DataFrame("", index=order, columns=tags)
    id_mat = pd.DataFrame("", index=order, columns=tags)
    flag_mat = np.zeros((n_rows, n_cols), dtype=bool)
    for i, fc in enumerate(order):
        present = top.get(fc, {})
        ids = [v[0] for v in present.values()]
        mode_id = Counter(ids).most_common(1)[0][0] if ids else None
        for j, tag in enumerate(tags):
            if tag in present:
                cid, r = present[tag]
                R_mat.loc[fc, tag] = r
                annot.loc[fc, tag] = f"{wrap_id(cid)}\n{r:.2f}"
                id_mat.loc[fc, tag] = cid
                flag_mat[i, j] = (cid != mode_id)

    core_h = ROW_PITCH_IN * n_rows
    core_w = ROW_PITCH_IN * n_cols
    fig_w = LEFT_MARGIN_IN + core_w + RIGHT_MARGIN_IN
    fig_h = TOP_MARGIN_IN + core_h + BOTTOM_MARGIN_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_MARGIN_IN / fig_w, BOTTOM_MARGIN_IN / fig_h, core_w / fig_w, core_h / fig_h])
    sns.heatmap(R_mat, ax=ax, cmap="Purples", vmin=0, vmax=np.nanmax(R_mat.values), square=True,
               linewidths=0.4, linecolor="white", cbar=False, xticklabels=True, yticklabels=True,
               annot=annot, fmt="", annot_kws={"fontsize": 6, "linespacing": 1.3},
               mask=R_mat.isna())
    ax.set_facecolor("#f2f2f2")
    cbar = fig.colorbar(ax.collections[0], ax=ax, shrink=0.45, aspect=12, pad=0.03)
    cbar.set_label("GSI Spearman r (top match)", fontsize=7)
    cbar.ax.tick_params(labelsize=LABEL_PT)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for tick, c in zip(ax.get_yticklabels(), R_mat.index):
        tick.set_fontsize(LABEL_PT)
        tick.set_color(row_color(c))
    ax.set_xticklabels(col_labels, fontsize=LABEL_PT, rotation=0, ha="center")
    ax.tick_params(axis="both", length=2, pad=1.5)

    for i in range(n_rows):
        for j in range(n_cols):
            if flag_mat[i, j]:
                ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor="#e31a1c", linewidth=1.6))

    present_pos = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order if c != SATB2_CLUSTER}]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present_pos]
    handles.append(Patch(facecolor=LILAC, edgecolor="none", label="SATB2-1"))
    handles.append(Patch(facecolor="none", edgecolor="#e31a1c", linewidth=1.6,
                         label="top match differs\nfrom mode across subsets"))
    heat_right_frac = (LEFT_MARGIN_IN + core_w) / fig_w
    heat_mid_frac = (BOTTOM_MARGIN_IN + core_h / 2) / fig_h
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(heat_right_frac + 0.16, heat_mid_frac),
              ncol=1, frameon=False, fontsize=6, title="row label colour", title_fontsize=6.5)

    fig.suptitle("GSI-SCT top chicken match per finch cluster, across subsets\n"
                "cell = matched chicken cluster + Spearman r; red border = identity differs from the cluster's own mode",
                fontsize=8)
    stub = "gsi_sct_subset_comparison"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220)
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")

    n_flagged = int(flag_mat.sum())
    print(f"\n{n_flagged} cell(s) flagged as composition-sensitive (identity != mode):")
    for i, fc in enumerate(order):
        for j, tag in enumerate(tags):
            if flag_mat[i, j]:
                print(f"  {fc} / {tag}: {id_mat.loc[fc, tag]} (r={R_mat.loc[fc, tag]:.3f})")


if __name__ == "__main__":
    main()
