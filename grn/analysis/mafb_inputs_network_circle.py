#!/usr/bin/env python3
"""Circle version of the MAFB-input network: the inputs to MAFB of BOTH cell groups (Glut-RA vs C1H-1 all runs, and GABA; config2) sit on one ellipse around
MAFB, RA inputs on the left half and GABA inputs on the right. Two panels share that layout: in the left panel the RA inputs are drawn bold and the GABA
inputs greyed out, in the right panel the reverse. Edge conventions as in mafb_inputs_network.py (arrow = positive TF2G rho, bar = negative, double head =
MAFB also targets the TF, dashed AR -> MAFB is hand-added: AR is not an eRegulon input to MAFB in the all-runs tables; self-links are not drawn).

    python mafb_inputs_network_circle.py [OUT_DIR]      # writes mafb_inputs_network_circle.{svg,pdf,png}
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

PANELS = [  # key, directory, title, background colour, input TFs along the arc, top to bottom (strongest TF-MAFB correlation first)
    ("RA_all", "glut-ra_glut-c1h-1_allruns_hybrid", "Glut-RA vs. Glut-C1H-1", "#e8504d", ["KLF12", "NR3C2", "NFIA", "KLF15", "LIN28B", "STAT1", "FOSL2"]),
    ("GABA", "ra-arco-hvc-nc_gaba_hybrid", "GABA", "#5b9a5b", ["MAF", "SOX6", "TCF12", "ARX", "SATB1", "MEIS2"]),
]
EXCLUDE = {}   # eRegulon inputs to leave out of the figure, {panel key: {TF}} (KLF15, weakly expressed in Glut-RA, is kept)
MANUAL = {"RA_all": ["AR"]}   # hand-added nodes (dashed edge into MAFB; see docstring)
INK, GREY, MAFB_FILL = "#1a1a1a", "#6b6b6b", "#d4d4d4"
cols = ["TF", "Gene", "regulation", "rho_TF2G"]


EXCLUDE = {}
MAFB_FILL = "#d4d4d4"
def edges_for(d):
    t = pd.concat([pd.read_table(f"{HB}{d}/results/config{CONFIG}/outs/{f}", usecols=cols) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
    inputs = set(t[(t.Gene == "MAFB") & (t.TF != "MAFB")].TF)
    nodes = inputs | {"MAFB"}
    e = t[t.TF.isin(nodes) & t.Gene.isin(nodes) & (t.TF != t.Gene)].groupby(["TF", "Gene"]).agg(rho=("rho_TF2G", "median"), reg=("regulation", "first")).reset_index()
    return inputs, e




import numpy as np
W, H = 7.0, 2.9                                # panel size (data units); lower H = more vertical compression
MX, MY = W / 2, H / 2                          # MAFB centre
RX, RY = 2.8, 1.2                              # ellipse half-axes
ANG = {"RA_all": (94, 266), "GABA": (68, -68)}   # first/last angle (degrees from +x, counter-clockwise) of each group's arc; RA left, GABA right
BG_ALPHA = 0.2                                   # opacity of the coloured panel backgrounds
STAGGER = 0.3                                  # every second node on a group's arc sits this fraction closer to MAFB (avoids overlaps in narrow panels)
GREY_FILL, GREY_INK, GREY_EDGE = "#c4c4c4", "#7a7a7a", "#9c9c9c"


def place(order, key):
    """nodes of one group on a half ellipse around MAFB, spaced evenly in y; every second node sits closer to MAFB (STAGGER) so neighbours do not collide"""
    a0, a1 = ANG[key]
    ys = np.linspace(np.sin(np.radians(a0)), np.sin(np.radians(a1)), len(order))
    sgn = np.sign(np.cos(np.radians(a0)))
    return {tf: (MX + sgn * RX * np.sqrt(1 - y * y) * (1 - STAGGER * (i % 2)), MY + RY * y) for i, (tf, y) in enumerate(zip(order, ys))}


def clip(c, hw, hh, toward, gap=0.05):
    dx, dy = toward[0] - c[0], toward[1] - c[1]
    s = min(hw / abs(dx) if dx else np.inf, hh / abs(dy) if dy else np.inf)
    n = np.hypot(dx, dy)
    return (c[0] + dx * s + gap * dx / n, c[1] + dy * s + gap * dy / n)


groups = {}                                    # key -> dict(order, nodes (incl. extras), edges)
for key, d, title, colour, order in PANELS:
    inputs, e = edges_for(d)
    inputs = inputs - EXCLUDE.get(key, set())
    assert set(order) == inputs, f"{key}: inputs changed: {inputs ^ set(order)}"
    nodes = order + MANUAL.get(key, [])
    groups[key] = dict(nodes=nodes, edges=e[((e.TF == "MAFB") | (e.Gene == "MAFB")) & (e.TF.isin(nodes + ["MAFB"])) & (e.Gene.isin(nodes + ["MAFB"]))],
                       title=title, colour=colour, pos=place(nodes, key), extra=MANUAL.get(key, []))
assert not (set(groups["RA_all"]["nodes"]) & set(groups["GABA"]["nodes"])), "the groups share an input TF: give the shared TF one position"
POS = {tf: p for g in groups.values() for tf, p in g["pos"].items()}; POS["MAFB"] = (MX, MY)
HH = 0.17
half = {tf: (0.165 * len(tf) + 0.3) / 2 for tf in POS}

fig, axes = plt.subplots(1, 2, figsize=(4.15, 1.5))
for ax, active in zip(axes, ("RA_all", "GABA")):
    A = groups[active]
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off"); ax.set_aspect("auto")
    ax.add_patch(FancyBboxPatch((0.15, 0.1), W - 0.3, H - 0.2, boxstyle="round,pad=0,rounding_size=0.45", facecolor=A["colour"], alpha=BG_ALPHA, edgecolor="none", zorder=0))
    ax.text(W / 2, H + 0.05, A["title"], ha="center", va="bottom", fontsize=6, fontweight="bold", clip_on=False)

    def ends(tf, shift=0.0):
        a, m = np.array(POS[tf]), np.array(POS["MAFB"])
        n = (m - a) / np.hypot(*(m - a)); off = shift * np.array([-n[1], n[0]])
        return clip(a + off, half[tf], HH, m + off), clip(m + off, half["MAFB"], HH, a + off)

    def arrow(start, end, head, on, dashed=False):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle=head, mutation_scale=1.0, linewidth=0.8 if on else 0.5, color=INK if on else GREY_EDGE,
                                     zorder=2 if on else 1, shrinkA=0, shrinkB=0, linestyle=(0, (3, 2)) if dashed else "-"))
    for key, g in groups.items():
        on = key == active
        pairs = set(zip(g["edges"].TF, g["edges"].Gene))
        for r in g["edges"].itertuples():
            if r.TF == "MAFB" and (r.Gene, "MAFB") in pairs:
                continue
            tf = r.TF if r.Gene == "MAFB" else r.Gene
            a_pt, m_pt = ends(tf)
            head = "-|>,head_width=2.2,head_length=3.2" if r.reg > 0 else "-[,widthB=3.2,lengthB=0.0"
            if r.Gene == "MAFB" and ("MAFB", tf) in pairs:
                head = "<|-|>,head_width=2.2,head_length=3.2"
            arrow(*((a_pt, m_pt) if r.Gene == "MAFB" else (m_pt, a_pt)), head, on)
        for tf in g["extra"]:                          # AR: dashed AR -> MAFB (hand-added), solid MAFB -> AR (config2 tables)
            a1, m1 = ends(tf, 0.07); a2, m2 = ends(tf, -0.07)
            arrow(a1, m1, "-|>,head_width=2.2,head_length=3.2", on, dashed=True)
            arrow(m2, a2, "-|>,head_width=2.2,head_length=3.2", on)
    for key, g in groups.items():
        for tf in g["nodes"]:
            x, y = POS[tf]; on = key == active; w = half[tf] * 2
            ax.add_patch(FancyBboxPatch((x - w / 2, y - HH), w, 2 * HH, boxstyle="round,pad=0,rounding_size=0.1", facecolor=g["colour"] if on else GREY_FILL,
                                        edgecolor=INK if on else GREY_EDGE, linewidth=0.8 if on else 0.5, zorder=3))
            ax.text(x, y, tf, ha="center", va="center", fontsize=6, style="italic", color=INK if on else GREY_INK, zorder=4)
    x, y = POS["MAFB"]; w = half["MAFB"] * 2
    ax.add_patch(FancyBboxPatch((x - w / 2, y - HH), w, 2 * HH, boxstyle="round,pad=0,rounding_size=0.1", facecolor=MAFB_FILL, edgecolor=INK, linewidth=0.8, zorder=3))
    ax.text(x, y, "MAFB", ha="center", va="center", fontsize=6, style="italic", color=INK, zorder=4)
fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.02, wspace=0.04)
for ext in ("svg", "pdf", "png"):
    fig.savefig(f"{OUT_DIR}/mafb_inputs_network_circle.{ext}", dpi=300)
print("wrote", OUT_DIR)
