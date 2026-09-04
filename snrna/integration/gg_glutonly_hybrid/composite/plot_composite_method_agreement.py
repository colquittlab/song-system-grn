"""Heatmap of per-method top-call agreement for the plusSATB2_noPre 3-method composite
(assemble_glutonly_plusSATB2_noPre.py's per_method_top_calls.csv) -- same fixed-inch
geometry/square-cell/row-colour/no-dendrogram conventions as the other heatmaps in this
analysis, adapted from plot_gsi_sct_subset_comparison.py's cell-annotation approach
(rows too few/text too long for plain colour-only cells).

Unlike the earlier GSI subset-comparison heatmap (continuous colour = correlation
magnitude), this one is CATEGORICAL: each cell's colour answers "does this method's own
top call for this finch cluster match the row's majority call?"
  - green  = matches the majority (either all 3 agree, or this is 2 of a 2-1 split)
  - orange = the odd one out in a 2-1 split (the one method that disagrees)
  - grey   = no majority at all (3-way disagreement -- "matches the mode" is meaningless
    here since every method picked something different)
Cell text: the matched chicken cluster's short name (Ex_ stripped, wrapped at the first
underscore to fit the cell).
"""
from pathlib import Path
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
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

GREEN = "#66c2a5"
ORANGE = "#fc8d62"
GREY = "#b3b3b3"

IN = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/composite/results/composite_plusSATB2_noPre")
OUT = IN

METHOD_COLS = ["gsi", "samap", "cca"]
METHOD_LABELS = {"gsi": "GSI", "samap": "SAMap", "cca": "CCA (RPCA)"}

BASE_ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]
ROW_ORDER = BASE_ROW_ORDER[:13] + [SATB2_CLUSTER] + BASE_ROW_ORDER[13:]

LABEL_PT = 6
ROW_PITCH_IN = 0.72
TOP_MARGIN_IN = 0.85
BOTTOM_MARGIN_IN = 0.55
LEFT_MARGIN_IN = 1.55
RIGHT_MARGIN_IN = 2.7


def row_color(c):
    return LILAC if c == SATB2_CLUSTER else POSITION_COLORS[CLUSTER_POSITION[c]]


def short(c):
    return c.replace("Ex_", "")


def wrap_id(cid):
    if "_" in cid and len(cid) > 8:
        fam, gene = cid.split("_", 1)
        return f"{fam}\n{gene}"
    return cid


def main():
    top = pd.read_csv(IN / "per_method_top_calls.csv", index_col=0)
    top = top.drop(index=[c for c in EXCLUDED_CLUSTERS if c in top.index])
    order = [c for c in ROW_ORDER if c in top.index]
    missing = set(top.index) - set(order)
    if missing:
        print(f"WARNING: clusters not in ROW_ORDER, appended at end: {sorted(missing)}")
        order += sorted(missing)
    top = top.loc[order]

    n_rows, n_cols = len(order), len(METHOD_COLS)
    code_mat = np.zeros((n_rows, n_cols), dtype=int)   # 0=green, 1=orange, 2=grey
    annot = pd.DataFrame("", index=order, columns=METHOD_COLS)
    for i, fc in enumerate(order):
        calls = top.loc[fc, METHOD_COLS]
        counts = Counter(calls)
        majority_n = counts.most_common(1)[0][1]
        for j, meth in enumerate(METHOD_COLS):
            cid = short(calls[meth])
            annot.loc[fc, meth] = wrap_id(cid)
            if majority_n == 1:
                code_mat[i, j] = 2  # total disagreement, no majority
            elif counts[calls[meth]] == majority_n:
                code_mat[i, j] = 0  # part of the majority
            else:
                code_mat[i, j] = 1  # odd one out

    cmap = ListedColormap([GREEN, ORANGE, GREY])
    core_h = ROW_PITCH_IN * n_rows
    core_w = ROW_PITCH_IN * n_cols
    fig_w = LEFT_MARGIN_IN + core_w + RIGHT_MARGIN_IN
    fig_h = TOP_MARGIN_IN + core_h + BOTTOM_MARGIN_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_MARGIN_IN / fig_w, BOTTOM_MARGIN_IN / fig_h, core_w / fig_w, core_h / fig_h])
    sns.heatmap(pd.DataFrame(code_mat, index=order, columns=METHOD_COLS), ax=ax, cmap=cmap,
               vmin=0, vmax=2, square=True, linewidths=0.4, linecolor="white", cbar=False,
               xticklabels=True, yticklabels=True, annot=annot, fmt="",
               annot_kws={"fontsize": 6, "linespacing": 1.3})
    for side in ax.spines.values():
        side.set_visible(True)
        side.set_edgecolor("black")
        side.set_linewidth(1.0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for tick, c in zip(ax.get_yticklabels(), order):
        tick.set_fontsize(LABEL_PT)
        tick.set_color(row_color(c))
    ax.set_xticklabels([METHOD_LABELS[m] for m in METHOD_COLS], fontsize=7, rotation=0, ha="center")
    ax.tick_params(axis="both", length=2, pad=1.5)

    present_pos = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in order if c != SATB2_CLUSTER}]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present_pos]
    if SATB2_CLUSTER in order:
        handles.append(Patch(facecolor=LILAC, edgecolor="none", label="SATB2-1"))
    agree_handles = [
        Patch(facecolor=GREEN, edgecolor="black", label="agrees w/ majority"),
        Patch(facecolor=ORANGE, edgecolor="black", label="odd one out (2-1 split)"),
        Patch(facecolor=GREY, edgecolor="black", label="no majority (3-way split)"),
    ]
    heat_right_frac = (LEFT_MARGIN_IN + core_w) / fig_w
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(heat_right_frac + 0.14, 0.85),
              ncol=1, frameon=False, fontsize=6.5, title="row label colour", title_fontsize=7)
    fig.legend(handles=agree_handles, loc="upper left", bbox_to_anchor=(heat_right_frac + 0.14, 0.45),
              ncol=1, frameon=False, fontsize=6.5, title="cell agreement", title_fontsize=7)

    n_full = int((code_mat == 0).all(axis=1).sum())
    n_partial = int((code_mat == 1).any(axis=1).sum())
    n_none = int((code_mat == 2).any(axis=1).sum())
    fig.suptitle("Per-method top-call agreement, plusSATB2_noPre 3-method composite (GSI, SAMap, CCA/RPCA)\n"
                f"{n_full}/{n_rows} full agreement, {n_partial}/{n_rows} 2-of-3, {n_none}/{n_rows} all disagree",
                fontsize=8)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / "method_agreement_heatmap.{}".format(ext), dpi=220)
    plt.close(fig)
    print(f"wrote {OUT}/method_agreement_heatmap.pdf / .png")


if __name__ == "__main__":
    main()
