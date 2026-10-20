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
import numpy as np
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
H = 4.2                                        # panel height (data units); lower H = more vertical compression
MX, MY = 8.3, H / 2                            # MAFB centre
RX, RY, SPAN = 5.8, 1.6, 68                   # the inputs sit on an ellipse around MAFB (half-axes, half-angle in degrees)


def arc(tfs):
    t = np.radians(np.linspace(-SPAN, SPAN, len(tfs)))
    return {tf: (MX - RX * np.cos(a), MY - RY * np.sin(a)) for tf, a in zip(tfs, t)}


def clip(c, hw, hh, toward, gap=0.05):
    """point where the line from box centre c to `toward` leaves the box (half-sizes hw, hh), plus a small gap"""
    dx, dy = toward[0] - c[0], toward[1] - c[1]
    s = min(hw / abs(dx) if dx else np.inf, hh / abs(dy) if dy else np.inf)
    n = np.hypot(dx, dy)
    return (c[0] + dx * s + gap * dx / n, c[1] + dy * s + gap * dy / n)



PANELS = [  # key, directory, title, background colour, input TFs along the arc, top to bottom (strongest TF-MAFB correlation first)
    ("RA_all", "glut-ra_glut-c1h-1_allruns_hybrid", "Glut-RA vs. Glut-C1H-1", "#e8504d", ["KLF12", "NR3C2", "NFIA", "KLF15", "LIN28B", "STAT1", "FOSL2"]),
    ("GABA", "ra-arco-hvc-nc_gaba_hybrid", "GABA", "#5b9a5b", ["MAF", "SOX6", "TCF12", "ARX", "SATB1", "MEIS2"]),
]
EXCLUDE = {}   # eRegulon inputs to leave out of the figure, {panel key: {TF}} (KLF15, weakly expressed in Glut-RA, is kept)
MANUAL = {"RA_all": ["AR"]}   # hand-added nodes (dashed edge into MAFB; see docstring)
INK, GREY, MAFB_FILL = "#1a1a1a", "#6b6b6b", "#d4d4d4"
cols = ["TF", "Gene", "regulation", "rho_TF2G"]


def edges_for(d):
    t = pd.concat([pd.read_table(f"{HB}{d}/results/config{CONFIG}/outs/{f}", usecols=cols) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
    inputs = set(t[(t.Gene == "MAFB") & (t.TF != "MAFB")].TF)
    nodes = inputs | {"MAFB"}
    e = t[t.TF.isin(nodes) & t.Gene.isin(nodes) & (t.TF != t.Gene)].groupby(["TF", "Gene"]).agg(rho=("rho_TF2G", "median"), reg=("regulation", "first")).reset_index()
    return inputs, e



fig, axes = plt.subplots(1, len(PANELS), figsize=(6.2, 2.0))
HH = 0.22                                      # box half-height
for ax, (key, d, title, colour, order) in zip(axes, PANELS):
    inputs, e = edges_for(d)
    inputs = inputs - EXCLUDE.get(key, set())
    assert set(order) == inputs, f"{key}: inputs changed: {inputs ^ set(order)}"
    extra = MANUAL.get(key, [])
    allnodes = order + extra
    pos = arc(allnodes); pos["MAFB"] = (MX, MY)
    e = e[((e.TF == "MAFB") | (e.Gene == "MAFB")) & e.TF.isin(pos) & e.Gene.isin(pos)]   # edges touching MAFB, between drawn nodes
    ax.set_xlim(0, 10); ax.set_ylim(0, H); ax.axis("off"); ax.set_aspect("auto")
    ax.add_patch(FancyBboxPatch((0.15, 0.1), 9.7, H - 0.2, boxstyle="round,pad=0,rounding_size=0.45", facecolor=colour, edgecolor="none", zorder=0))
    ax.text(5, H + 0.05, title, ha="center", va="bottom", fontsize=6, fontweight="bold", clip_on=False)
    half = {tf: (0.15 * len(tf) + 0.4) / 2 for tf in pos}
    pairs = set(zip(e.TF, e.Gene))

    def ends(tf, shift=0.0):
        """(point on the TF box, point on the MAFB box) for the straight edge between them; shift offsets the line sideways (data units)"""
        a, m = np.array(pos[tf]), np.array(pos["MAFB"])
        n = (m - a) / np.hypot(*(m - a)); off = shift * np.array([-n[1], n[0]])
        return clip(a + off, half[tf], HH, m + off), clip(m + off, half["MAFB"], HH, a + off)

    def arrow(start, end, head, dashed=False):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle=head, mutation_scale=1.0, linewidth=0.8, color=INK, zorder=2, shrinkA=0, shrinkB=0,
                                     linestyle=(0, (3, 2)) if dashed else "-"))
    for r in e.itertuples():
        if r.TF == "MAFB" and (r.Gene, "MAFB") in pairs:
            continue                                  # drawn once with the TF -> MAFB edge, as a double head
        tf = r.TF if r.Gene == "MAFB" else r.Gene
        a_pt, m_pt = ends(tf)
        head = "-|>,head_width=2.2,head_length=3.2" if r.reg > 0 else "-[,widthB=3.2,lengthB=0.0"
        if r.Gene == "MAFB" and ("MAFB", tf) in pairs:
            head = "<|-|>,head_width=2.2,head_length=3.2"
        start, end = (a_pt, m_pt) if r.Gene == "MAFB" else (m_pt, a_pt)   # MAFB -> TF only: reverse the direction
        arrow(start, end, head)
    for tf in extra:                                   # AR: dashed AR -> MAFB (not in the tables), solid MAFB -> AR (config2 tables)
        a1, m1 = ends(tf, 0.07); a2, m2 = ends(tf, -0.07)
        arrow(a1, m1, "-|>,head_width=2.2,head_length=3.2", dashed=True)
        arrow(m2, a2, "-|>,head_width=2.2,head_length=3.2")
    if extra:
        ax.text(0.5, 0.22, "dashed: not an eRegulon edge in this run", ha="left", va="bottom", fontsize=4.5, color=INK)
    for tf, (x, y) in pos.items():
        w = half[tf] * 2
        ax.add_patch(FancyBboxPatch((x - w / 2, y - HH), w, 2 * HH, boxstyle="round,pad=0,rounding_size=0.14",
                                    facecolor=MAFB_FILL if tf == "MAFB" else colour, edgecolor=INK, linewidth=0.8, zorder=3))
        ax.text(x, y, tf, ha="center", va="center", fontsize=6, style="italic", color=INK, zorder=4)
fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.02, wspace=0.04)
for ext in ("svg", "pdf", "png"):
    fig.savefig(f"{OUT_DIR}/mafb_inputs_network.{ext}", dpi=300)
print("wrote", OUT_DIR)
