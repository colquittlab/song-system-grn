#!/usr/bin/env python3
"""Network mock-up of the TFs that regulate MAFB (config2) in RA vs C1H-1 (all runs) and in GABA, in the style of the earlier Figure 6 panel
(rounded coloured panel per cell group, italic TF boxes, MAFB shaded grey).

Nodes: MAFB plus every TF that is an input to MAFB in the dataset's config2 eRegulon tables (direct + extended). Edges: only those touching MAFB (TF -> MAFB, and\nMAFB -> TF where the same tables have it; a double head marks both). Arrow = positive regulation (TF2G rho > 0), flat bar = negative.
AR (RA panel) is added by hand: it is not an eRegulon input to MAFB in the all-runs config1/2 tables (its direct cistrome equals NR3C2's), so AR -> MAFB is
drawn dashed; MAFB -> AR is in the config2 tables (19 rows, positive) and is drawn solid.
TF self-links (TF2G rho exactly 1.000) are injected by SCENIC+ and are not drawn.

    python mafb_inputs_network.py [OUT_DIR]      # writes mafb_inputs_network.{svg,pdf,png}; SVG text stays editable
"""
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

warnings.filterwarnings("ignore")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

HB = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/ssd/rstudio/multiome/motor-pathway/scenicplus/mafb_inputs")
os.makedirs(OUT_DIR, exist_ok=True)
CONFIG = 2
def column(tfs, x=2.4, top=5.5, bottom=0.62):
    step = (top - bottom) / (len(tfs) - 1)
    return {tf: (x, top - i * step) for i, tf in enumerate(tfs)}


PANELS = [  # key, directory, title, background colour, input TFs top to bottom (strongest TF-MAFB correlation first)
    ("RA_all", "glut-ra_glut-c1h-1_allruns_hybrid", "Glut-RA vs. Glut-C1H-1", "#e8504d", ["KLF12", "NR3C2", "NFIA", "KLF15", "LIN28B", "STAT1", "FOSL2"]),
    ("GABA", "ra-arco-hvc-nc_gaba_hybrid", "GABA", "#5b9a5b", ["MAF", "SOX6", "TCF12", "ARX", "SATB1", "MEIS2"]),
]
MANUAL = {"RA_all": ["AR"]}   # hand-added nodes (dashed edge into MAFB; see docstring)
INK, GREY, MAFB_FILL = "#1a1a1a", "#6b6b6b", "#d4d4d4"
cols = ["TF", "Gene", "regulation", "rho_TF2G"]


def edges_for(d):
    t = pd.concat([pd.read_table(f"{HB}{d}/results/config{CONFIG}/outs/{f}", usecols=cols) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
    inputs = set(t[(t.Gene == "MAFB") & (t.TF != "MAFB")].TF)
    nodes = inputs | {"MAFB"}
    e = t[t.TF.isin(nodes) & t.Gene.isin(nodes) & (t.TF != t.Gene)].groupby(["TF", "Gene"]).agg(rho=("rho_TF2G", "median"), reg=("regulation", "first")).reset_index()
    return inputs, e


fig, axes = plt.subplots(1, len(PANELS), figsize=(6.2, 2.5))
for ax, (key, d, title, colour, order) in zip(axes, PANELS):
    inputs, e = edges_for(d)
    assert set(order) == inputs, f"{key}: inputs changed: {inputs ^ set(order)}"
    extra = MANUAL.get(key, [])
    allnodes = order + extra
    pos = column(allnodes); pos["MAFB"] = (8.1, 3.0)
    e = e[(e.TF == "MAFB") | (e.Gene == "MAFB")]
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off"); ax.set_aspect("auto")
    ax.add_patch(FancyBboxPatch((0.15, 0.1), 9.7, 5.8, boxstyle="round,pad=0,rounding_size=0.45", facecolor=colour, edgecolor="none", zorder=0))
    ax.text(5, 6.15, title, ha="center", va="bottom", fontsize=6, fontweight="bold", clip_on=False)
    half = {tf: (0.17 * len(tf) + 0.45) / 2 for tf in pos}
    pairs = set(zip(e.TF, e.Gene))
    mx, my = pos["MAFB"]
    slot = {tf: my + 0.62 * (i / (len(allnodes) - 1) - 0.5) * -1 for i, tf in enumerate(allnodes)}   # spread the arrow ends along MAFB's left edge
    for r in e.itertuples():
        if r.TF == "MAFB" and (r.Gene, "MAFB") in pairs:
            continue                                  # drawn once with the TF -> MAFB edge, as a double head
        tf = r.TF if r.Gene == "MAFB" else r.Gene
        a_pt = (pos[tf][0] + half[tf] + 0.04, pos[tf][1]); m_pt = (mx - half["MAFB"] - 0.04, slot[tf])
        both = r.Gene == "MAFB" and ("MAFB", tf) in pairs
        head = "-|>,head_width=2.2,head_length=3.2" if r.reg > 0 else "-[,widthB=3.2,lengthB=0.0"
        if both:
            head = "<|-|>,head_width=2.2,head_length=3.2"
        start, end = (a_pt, m_pt) if r.Gene == "MAFB" else (m_pt, a_pt)   # MAFB -> TF only: reverse the direction
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle=head, mutation_scale=1.0, linewidth=0.8, color=INK, zorder=2, shrinkA=0, shrinkB=0))
    for tf in extra:                                   # AR: dashed AR -> MAFB (not in the tables), solid MAFB -> AR (config2 tables)
        a_pt = (pos[tf][0] + half[tf] + 0.04, pos[tf][1]); m_pt = (mx - half["MAFB"] - 0.04, slot[tf])
        up, down = (m_pt[0], m_pt[1] + 0.09), (m_pt[0], m_pt[1] - 0.09)
        ax.add_patch(FancyArrowPatch((a_pt[0], a_pt[1] + 0.05), up, arrowstyle="-|>,head_width=2.2,head_length=3.2", mutation_scale=1.0, linewidth=0.8,
                                     color=INK, linestyle=(0, (3, 2)), zorder=2, shrinkA=0, shrinkB=0))
        ax.add_patch(FancyArrowPatch(down, (a_pt[0], a_pt[1] - 0.05), arrowstyle="-|>,head_width=2.2,head_length=3.2", mutation_scale=1.0, linewidth=0.8,
                                     color=INK, zorder=2, shrinkA=0, shrinkB=0))
    if extra:
        ax.text(9.7, 0.3, "dashed: not an eRegulon edge in this run", ha="right", va="bottom", fontsize=4.5, color=INK)
    for tf, (x, y) in pos.items():
        w = half[tf] * 2
        ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.31), w, 0.62, boxstyle="round,pad=0,rounding_size=0.2",
                                    facecolor=MAFB_FILL if tf == "MAFB" else colour, edgecolor=INK, linewidth=0.8, zorder=3))
        ax.text(x, y, tf, ha="center", va="center", fontsize=6, style="italic", color=INK, zorder=4)
fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.02, wspace=0.04)
for ext in ("svg", "pdf", "png"):
    fig.savefig(f"{OUT_DIR}/mafb_inputs_network.{ext}", dpi=300)
print("wrote", OUT_DIR)
