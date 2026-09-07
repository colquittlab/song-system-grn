"""Hierarchically-clustered heatmap of a composite rank_score matrix, class-annotated.

Design decisions and why:

MATRIX COLOUR — rank_score is a magnitude in [0, 1] (1 = every method ranked that
reference label first), so the ramp is SEQUENTIAL and single-hue, light->dark. Orange is
used deliberately because it is NOT one of the six categorical hues in the annotation
strips, so a dark cell can never be confused with a class colour.

CLASS COLOURS — validated with the dataviz palette validator under --pairs all (any two
classes may end up compared, not just neighbours). The documented 8-slot categorical
palette FAILS all-pairs at 8, 6, 5 and 4 slots on the normal-vision floor; a search over
subsets found 6 to be the maximum that passes, hence the eight biological classes are
merged to six (astro/oligo/ependymal -> "glia"). Chosen set passes with:
    normal-vision floor 15.6 (>=15 required), CVD worst 6.9 (6-8 band).
The CVD 6-8 band is legal ONLY with secondary encoding, and the contrast WARN on aqua and
yellow requires visible relief, so BOTH are provided: a legend, and a bracketed class tag
appended to every tick label ([neu], [nbl], [pro], [gli], [imm], [vas]) so class is always
readable as text and never colour-alone.

COLUMN FILTERING — a reference has 400-550 labels, far more than can be rendered
legibly. Columns are reduced to the union of each finch cluster's top-K matches. The
number dropped is printed and written into the figure subtitle rather than silently
truncated, since "we plotted everything" and "we plotted the top slice" look identical
on the page otherwise.
"""
import argparse, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

# Established project convention (analyze_integration.py): signed data (a correlation,
# not a magnitude) gets a diverging ramp through white, never a sequential one -- a
# sequential ramp on signed data makes a strong negative and a weak positive look alike.
DIVERGING = LinearSegmentedColormap.from_list("corr_div", ["#2a78d6", "#ffffff", "#eb6834"])
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from class_benchmark import expected_class, ref_class

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})


def dot_offsets(n: int):
    """(dx, dy) offsets, in units of one cell (cell spans [-0.5, 0.5] on each axis), for
    n agreement-count dots -- Zaremba et al.'s convention of one dot per method agreeing
    on a pair. n<=3 lays out a single centred row (matches the visual style of the
    published figure); n>=4 falls back to a roughly square grid so dots stay legible
    inside a small square cell instead of over-stretching a row."""
    import math
    if n <= 0:
        return []
    if n == 1:
        return [(0.0, 0.0)]
    if n == 2:
        span = 0.16                          # pulled in tighter toward centre
        return [(-span, 0.0), (span, 0.0)]
    if n == 3:
        span = 0.30
        return [(x, 0.0) for x in np.linspace(-span, span, 3)]
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    span = 0.20                              # pulled in tighter toward centre
    xs = np.linspace(-span, span, ncols)
    ys = np.linspace(-span, span, nrows)
    return [(x, y) for y in ys for x in xs][:n]

LABEL_PT = 6.5   # tick-label size; cell pitch is derived from it

# 6-class merge of the 8 benchmark classes (see module docstring), EXCEPT "neuron" is
# further split into "glut"/"gaba" (excitatory/inhibitory) by a name-pattern check
# (glut_gaba_split, below) rather than by expected_class/ref_class -- those two stay
# generic ("neuron") because class_benchmark.py's accuracy scoring depends on that coarser
# grouping and must not change; this split is display-only, local to this script.
MERGE = {"neuron": "neuron", "neuroblast": "neuroblast", "progenitor": "progenitor",
         "astro": "glia", "oligo": "glia", "ependymal": "glia",
         "immune": "immune", "vascular": "vascular", "unknown": "unknown"}
ORDER = ["glut", "gaba", "neuroblast", "progenitor", "glia", "immune", "vascular"]
# validated all-pairs, light surface (node scripts/validate_palette.js, dataviz skill):
# normal-vision floor 15.6, worst CVD (deutan) 6.9 -- in the 6-8 floor band, legal only
# with secondary encoding, provided here by the legend + the strip position itself.
COLORS = {"glut": "#2a78d6", "gaba": "#a34e9e", "neuroblast": "#1baf7a", "progenitor": "#eda100",
          "glia": "#008300", "immune": "#4a3aa7", "vascular": "#e34948",
          "unknown": "#c9c9c4", "neuron": "#2a78d6"}   # neuron: safety-net alias of glut,
          # for a name classified "neuron" that glut_gaba_split can't further resolve


def glut_gaba_split(name: str, fallback: str) -> str:
    """Excitatory/inhibitory split by name pattern -- covers every vocabulary seen in this
    project so far: finch/mouse 'Glut-'/'GABA-' prefixes, chicken 'Ex_'/'Inh_' prefixes.
    Anything not clearly one or the other keeps whatever coarse class it already had
    (only meaningful for the "neuron" bucket; glia/immune/etc. are untouched)."""
    n = name.lower()
    if n.startswith("gaba") or n.startswith("inh_") or n.startswith("inh-") or "gaba" in n:
        return "gaba"
    if n.startswith("glut") or n.startswith("ex_") or n.startswith("ex-") or "glut" in n:
        return "glut"
    return fallback

# Validated all-pairs, light surface: CVD 21.6, normal-vision 32.3. Used in place of the
# coarse-class strip when every row is the same class (e.g. a Glut-only subset), where a
# monochrome class strip would carry no information.
GROUP_COLORS = {"song": "#e34948", "non-song": "#2a78d6"}


def main(matrix: Path, out_prefix: Path, annot_csv: Path | None, top_k: int,
         title: str, max_cols: int, cbar_label: str,
         cell: float | None, group_csv: Path | None, signed: bool = False,
         dots_csv: Path | None = None, dots_min: int = 2,
         transpose: bool = False, scale: float = 1.0, no_dendrograms: bool = False,
         label_pt_override: float | None = None, row_order_csv: Path | None = None):
    M = pd.read_csv(matrix, index_col=0)
    annot = pd.read_csv(annot_csv, index_col=0) if annot_csv and annot_csv.exists() else None
    n_all = M.shape[1]

    # load the full agreement-count matrix (if given) BEFORE column selection, so each
    # row's peak-agreement column can be preserved even when it isn't in that row's
    # top-K by confidence -- see the GABA-5-2 case: top-1 confidence beat the 3-method
    # reciprocal-consensus column by only 0.005, which top_k=1 alone would silently drop.
    dots_full = None
    if dots_csv and dots_csv.exists():
        dd = pd.read_csv(dots_csv, index_col=0)
        dd.columns = dd.columns.astype(str)
        dots_full = dd.reindex(index=M.index, columns=M.columns).fillna(0).astype(int)

    # keep the union of each row's top-K columns by confidence...
    keep = set()
    for c in M.index:
        keep |= set(M.loc[c].nlargest(top_k).index)
    # ...plus each row's own best-agreement column, so a column is never dropped while
    # still holding that row's peak cross-method consensus (only meaningful at >=dots_min,
    # matching the threshold dots are actually drawn at).
    anchor_cols = set()
    if dots_full is not None:
        for c in M.index:
            ra = dots_full.loc[c]
            if ra.max() >= dots_min:
                anchor_cols.add(ra.idxmax())
    keep |= anchor_cols
    keep = [c for c in M.columns if c in keep]
    if len(keep) > max_cols:                       # fall back to global strength, but never
                                                    # trim away an agreement anchor column
        trimmable = [c for c in keep if c not in anchor_cols]
        n_trimmable_keep = max(max_cols - len(anchor_cols), 0)
        strength = M[trimmable].max(axis=0).nlargest(n_trimmable_keep) if trimmable else pd.Series(dtype=float)
        keep = [c for c in M.columns if c in anchor_cols or c in set(strength.index)]
    Msub = M[keep]
    dropped = n_all - len(keep)
    anchor_note = f"; +{len(anchor_cols)} agreement-anchor cols (>= {dots_min})" if anchor_cols else ""
    print(f"{matrix.name}: {M.shape[0]} x {n_all} -> kept {len(keep)} columns "
          f"(union of per-cluster top-{top_k}{anchor_note}); dropped {dropped}")

    # optional FIXED row order (e.g. matching the biological ordering used by this
    # project's other, non-clustered heatmaps), replacing hierarchical row clustering.
    # Rows in Msub but not in the file are DROPPED (not appended) -- the file is the
    # authoritative row set, so this also serves as an exclusion list (e.g. an edge-case
    # cluster this analysis excludes everywhere else).
    row_order = None
    if row_order_csv and row_order_csv.exists():
        row_order = pd.read_csv(row_order_csv, header=None)[0].tolist()
        keep_rows = [c for c in row_order if c in Msub.index]
        dropped_rows = set(Msub.index) - set(keep_rows)
        if dropped_rows:
            print(f"  row_order_csv: dropping {len(dropped_rows)} row(s) not listed: {sorted(dropped_rows)}")
        Msub = Msub.loc[keep_rows]

    # optional song/non-song (or any other) row grouping, REPLACING the coarse-class
    # strip -- meant for subsets (e.g. Glut-only) where every row shares one class and a
    # class strip would be uninformative.
    group = None
    if group_csv and group_csv.exists():
        gg = pd.read_csv(group_csv, index_col=0)
        if "group" in gg.columns:
            group = gg["group"].reindex(Msub.index)

    # optional per-cell cross-method agreement count, overlaid as dots (Zaremba
    # convention): N dots in a cell = N methods independently call that pair a
    # reciprocal top-N match. Aligned now, on Msub's (pre-clustering) row/col order;
    # permuted into the post-clustering plotted order once the dendrograms exist below.
    dots = dots_full.reindex(index=Msub.index, columns=Msub.columns) if dots_full is not None else None

    if group is not None:
        row_cls = list(group.fillna("unknown"))
        ROW_COLORS_MAP = {**GROUP_COLORS, "unknown": "#c9c9c4"}
    else:
        row_cls = [MERGE.get(expected_class(c) or "unknown", "unknown") for c in Msub.index]
        row_cls = [glut_gaba_split(c, cls) if cls == "neuron" else cls
                  for c, cls in zip(Msub.index, row_cls)]
        ROW_COLORS_MAP = COLORS
    col_cls = [MERGE.get(ref_class(c, annot), "unknown") for c in Msub.columns]  # column strip stays coarse-class regardless
    col_cls = [glut_gaba_split(c, cls) if cls == "neuron" else cls
              for c, cls in zip(Msub.columns, col_cls)]
    row_colors = pd.DataFrame({"class": [ROW_COLORS_MAP[c] for c in row_cls]}, index=Msub.index)
    col_colors = pd.Series([COLORS[c] for c in col_cls], index=Msub.columns, name="class")

    if transpose:
        # Flip which axis displays which -- done AFTER row_cls/col_cls/row_colors/
        # col_colors are computed with their correct semantic meaning (expected_class for
        # finch, ref_class for the reference species), so this is a pure presentational
        # swap, not a re-classification.
        Msub = Msub.T
        row_colors, col_colors = (pd.DataFrame({"class": col_colors.values}, index=col_colors.index),
                                  pd.Series(row_colors["class"].values, index=row_colors.index, name="class"))
        if dots is not None:
            dots = dots.T

    # Class/group identity is carried by the colour strips + legend only. (An earlier
    # version also appended a bracketed tag to every tick label as a secondary,
    # colour-independent encoding for CVD accessibility -- removed per explicit request;
    # the strips + legend remain the source of truth for class/confidence.)
    Mp = Msub.copy()
    row_colors.index = Mp.index
    col_colors.index = Mp.columns

    # Cell size is set from the label font size so tick labels sit adjacent with minimal
    # whitespace: at LABEL_PT the text occupies LABEL_PT/72 in, so a cell only slightly
    # larger than that packs labels tightly without overlapping. A single SQUARE pitch is
    # used for both axes -- pitch_row/pitch_col equal the actual rendered cell height/width
    # in inches (see h_core/w below), so setting them equal makes every cell square
    # regardless of row/column count.
    # --scale shrinks (or grows) the WHOLE figure proportionally -- every absolute size
    # below (font points, strip/margin inches, dot diameter) is multiplied by it, so a
    # scale=0.2 render is a true miniature of the full design, not the same fonts crammed
    # into a smaller canvas (which is what plainly resizing the saved figure would give).
    # --label_pt_override decouples the RENDERED TICK LABEL size from that scaling -- e.g.
    # scale=0.8 (a 20% smaller figure) with label_pt_override=6 still derives cell/strip/
    # margin geometry from the scaled size (so the figure is genuinely 20% smaller), but
    # renders tick labels at a fixed, legible point size instead of shrinking them too.
    label_pt = LABEL_PT * scale
    render_label_pt = label_pt_override if label_pt_override else label_pt
    base = (label_pt / 72.0)
    pitch_row = cell if cell else base * 1.65   # extra headroom: real glyph line-height (ascenders/descenders, "_") exceeds nominal font size
    pitch_col = pitch_row
    # The pitch must apply to the HEATMAP AXES, not the whole figure: clustermap lays the
    # axes out as a FRACTION of the figure (the dendrograms and colour strips take the
    # rest), so sizing the figure directly by pitch*n squeezed each row below 5pt and the
    # labels collided. Divide by the axes fraction to get the figure size.
    # Hierarchical clustering still runs (and still sets row/column order) even with the
    # dendrogram lines hidden -- only the drawn tree and its reserved gridspec space go
    # away, freeing that room for the heatmap itself. A ratio of exactly 0 isn't accepted
    # by seaborn, so use a negligible sliver instead.
    dr_row, dr_col = (0.001, 0.001) if no_dendrograms else (0.10, 0.06)  # dendrogram_ratio below
    if row_order is not None:
        dr_row = 0.001  # no row clustering happens at all -- nothing to reserve tree space for
    # STRIP_IN is the class-colour strip's thickness in ABSOLUTE inches, held equal on
    # both axes so the left (row) and top (column) panels read as the same visual weight
    # regardless of the matrix's aspect ratio. colors_ratio is a FRACTION of each axis's
    # own total (w for the row strip's width, h_core for the column strip's height), so
    # a single numeric ratio does not give equal absolute thickness when w != h_core --
    # solved in closed form instead: w = pitch*n_cols/(1-dr_row-cr_row) with
    # cr_row = STRIP_IN/w rearranges to w = (pitch*n_cols + STRIP_IN)/(1-dr_row), and
    # symmetrically for h_core/cr_col.
    STRIP_IN = 0.09 * scale
    # Suptitle + legend sit above the clustermap axes via bbox/y>1 anchors, at y-fractions
    # of the TOTAL figure height. Those fractions were tuned for tall (11x5in+) figures; on
    # a short one (few rows, e.g. a 20-row subset) the same fractions compress into too
    # little ABSOLUTE space and the title/legend/column-dendrogram collide. Fix: add a
    # fixed-INCH margin band on top of h_core (not a fraction of it), push the clustermap's
    # own axes down into the bottom frac_top of the enlarged canvas via subplots_adjust, and
    # anchor title/legend within the newly-freed band -- so every element gets the same
    # absolute room regardless of row count.
    MARGIN_IN = 1.25 * scale   # bumped: legend can now wrap to 2 rows (ncol capped at 4)
    cmap = DIVERGING if signed else "Oranges"
    vmin, vmax = (-1, 1) if signed else (0, 1)
    CBAR_POS = (0.91, 0.02, 0.006, 0.11)   # lower-right, narrow; figure-fraction (left, bottom, width, height)
    n_rows, n_cols = Mp.shape[0], Mp.shape[1]

    # seaborn's actual rendered ax_heatmap size does not match the analytical
    # dendrogram_ratio/colors_ratio/subplots_adjust math closely enough to trust blindly
    # (measured to be off by ~15-35% in both axes on real matrices) -- so MEASURE the
    # real cell size after building, then rebuild once with a corrected figure size. Two
    # passes converge to <1% because the relationship is close to linear for a small
    # correction.
    w = max(5.0 * scale, pitch_col * n_cols / (1.0 - dr_row - 0.012))
    h_core = max(3.0 * scale, pitch_row * n_rows / (1.0 - dr_col - 0.012))
    g = None
    for _pass in range(8):
        cr_row = STRIP_IN / w
        cr_col = STRIP_IN / h_core
        h = h_core + MARGIN_IN
        frac_top = h_core / h
        if g is not None:
            plt.close(g.figure)
        g = sns.clustermap(
            Mp, cmap=cmap, vmin=vmin, vmax=vmax,
            row_colors=row_colors, col_colors=col_colors,
            figsize=(w, h), linewidths=0, rasterized=True,
            dendrogram_ratio=(dr_row, dr_col), colors_ratio=(cr_row, cr_col),
            cbar_pos=CBAR_POS,
            xticklabels=True, yticklabels=True,
            row_cluster=(row_order is None),
            metric="euclidean", method="average",
        )
        g.figure.subplots_adjust(top=frac_top)
        bbox = g.ax_heatmap.get_position()
        fw, fh = g.figure.get_size_inches()
        cell_w = bbox.width * fw / n_cols
        cell_h = bbox.height * fh / n_rows
        if abs(cell_w - pitch_col) < 0.001 and abs(cell_h - pitch_row) < 0.001:
            break
        w *= pitch_col / cell_w
        h_core *= pitch_row / cell_h
    # subplots_adjust resets any axes still tracked by the figure's gridspec back to its
    # own computed position -- including ax_cbar, silently undoing cbar_pos above. Re-apply
    # it now that subplots_adjust has already run.
    g.ax_cbar.set_position(CBAR_POS)
    if no_dendrograms and g.ax_row_dendrogram is not None:
        g.ax_row_dendrogram.set_visible(False)
    if no_dendrograms and g.ax_col_dendrogram is not None:
        g.ax_col_dendrogram.set_visible(False)
    # seaborn labels each colour strip with its row_colors/col_colors column name
    # ("class") by default -- redundant with the legend, drop it.
    if g.ax_row_colors is not None:
        g.ax_row_colors.set_xticklabels([]); g.ax_row_colors.set_xlabel("")
    if g.ax_col_colors is not None:
        g.ax_col_colors.set_yticklabels([]); g.ax_col_colors.set_ylabel("")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90,
                                 fontsize=render_label_pt)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0,
                                 fontsize=render_label_pt)
    # trim the gap between the axes and the tick text
    g.ax_heatmap.tick_params(axis="both", which="both", length=1.0 * scale, pad=0.8 * scale)
    g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("")
    # label ABOVE the bar: as a rotated ylabel it overlapped its own tick numbers
    g.cax.set_ylabel("")
    g.cax.set_title(cbar_label, fontsize=6.5 * scale, pad=3 * scale)
    g.cax.set_yticks([vmin, vmax])   # only the scale's min/max, not intermediate ticks
    g.cax.tick_params(labelsize=6 * scale)

    if dots is not None:
        # clustermap reorders both axes by its own hierarchical clustering; dots must be
        # placed at the PLOTTED position, not Msub's original row/column order. With
        # row_cluster=False (row_order given), g.dendrogram_row is None -- rows are
        # already in the plotted (fixed) order, so no permutation is needed there.
        row_ind = g.dendrogram_row.reordered_ind if g.dendrogram_row is not None else np.arange(Mp.shape[0])
        col_ind = g.dendrogram_col.reordered_ind if g.dendrogram_col is not None else np.arange(Mp.shape[1])
        Dp = dots.iloc[row_ind, col_ind]
        dot_diam_pt = pitch_row * 72.0 * 0.20   # smaller: 0.34 overlapped neighbouring cells
        xs, ys = [], []
        for ridx in range(Dp.shape[0]):
            for cidx in range(Dp.shape[1]):
                n = int(Dp.iat[ridx, cidx])
                if n < dots_min:
                    continue
                for dx, dy in dot_offsets(n):
                    xs.append(cidx + 0.5 + dx)
                    ys.append(ridx + 0.5 + dy)
        if xs:
            g.ax_heatmap.scatter(xs, ys, s=dot_diam_pt ** 2, c="black",
                                 edgecolors="none", zorder=10)

    if group is not None:
        row_present = [c for c in ["song", "non-song", "unknown"] if c in set(row_cls)]
        handles = [Patch(facecolor=ROW_COLORS_MAP[c], edgecolor="none", label=c) for c in row_present]
        col_present = [c for c in ORDER if c in set(col_cls)]
        handles += [Patch(facecolor=COLORS[c], edgecolor="none", label=f"{c} (chicken class)")
                   for c in col_present]
    else:
        present = [c for c in ORDER if c in set(row_cls) | set(col_cls)]
        handles = [Patch(facecolor=COLORS[c], edgecolor="none", label=c) for c in present]
        if "unknown" in set(row_cls) | set(col_cls):
            handles.append(Patch(facecolor=COLORS["unknown"], edgecolor="none", label="unknown"))
    # Bottom-centre, figure-level: anchoring to ax_col_dendrogram put the legend on top of
    # seaborn's own "reference class" strip label at the right edge.
    g.figure.legend(handles=handles, loc="lower center",
                    bbox_to_anchor=(0.5, frac_top + 0.10 * (1 - frac_top)),
                    # A single-row legend (ncol=len(handles)) can render far wider than the
                    # nominal figure width when there are many handles; bbox_inches="tight"
                    # then expands the canvas to fit it, which distorts the aspect ratio and
                    # squashes row pitch. Capping at 4 columns keeps the legend's rendered
                    # width in the same ballpark as the heatmap regardless of handle count.
                    frameon=False, fontsize=7 * scale, ncol=min(len(handles), 4),
                    title=("row/column grouping (colour strips)" if group is not None
                          else "coarse class (colour strips)"),
                    title_fontsize=7.5 * scale)
    # State the observed range explicitly: the colour scale is fixed at 0-1 so panels are
    # comparable ACROSS comparisons, which means a genuinely weak comparison looks pale.
    # Without this note a pale panel is indistinguishable from a truncated colour scale.
    scale_note = "-1 to 1 (diverging, signed)" if signed else "0-1"
    dot_note = (f"   ·   dots = methods agreeing this pair is a reciprocal top-N match "
               f"(Zaremba); shown for ≥{dots_min}" if dots is not None else "")
    sub = (f"{cbar_label.replace(chr(10), ' ')} — colour scale fixed {scale_note} for cross-panel "
           f"comparability; observed range here {Msub.values.min():.2f}–{Msub.values.max():.2f}"
           f"   ·   {len(keep)} of {n_all} reference labels shown "
           f"(union of per-cluster top-{top_k}; {dropped} omitted){dot_note}")
    g.figure.suptitle(f"{title}\n{sub}", fontsize=7.5 * scale, x=0.5,
                      y=frac_top + 0.85 * (1 - frac_top), ha="center")
    for ext in ("pdf", "png"):
        g.figure.savefig(f"{out_prefix}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(g.figure)
    print(f"  wrote {out_prefix}.pdf / .png  ({w:.0f}x{h:.0f} in)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--out_prefix", type=Path, required=True)
    p.add_argument("--annot_csv", type=Path, default=None)
    p.add_argument("--top_k", type=int, default=3)
    p.add_argument("--max_cols", type=int, default=90)
    p.add_argument("--title", default="composite rank score")
    p.add_argument("--cbar_label", default="composite\nrank score")
    p.add_argument("--group_csv", type=Path, default=None,
                   help="CSV with a 'group' column indexed by finch cluster (e.g. song/non-song); "
                        "replaces the coarse-class row strip when given")
    p.add_argument("--signed", action="store_true",
                   help="Matrix is a signed correlation (e.g. GSI), not a [0,1] magnitude: "
                        "use a diverging blue-white-orange ramp over [-1,1] instead of sequential Oranges over [0,1]")
    p.add_argument("--cell", type=float, default=None,
                   help="inches per row/column; default derives from the 5pt label size")
    p.add_argument("--dots_csv", type=Path, default=None,
                   help="composite_agreement_count_matrix.csv; overlays N dots per cell "
                        "for N methods agreeing (Zaremba convention)")
    p.add_argument("--dots_min", type=int, default=2,
                   help="minimum agreement count to draw dots for (default 2: 1 method isn't 'agreement')")
    p.add_argument("--transpose", action="store_true",
                   help="swap which axis displays rows vs columns (presentational only -- "
                        "classification/filtering is computed before the flip)")
    p.add_argument("--scale", type=float, default=1.0,
                   help="uniformly scale the whole figure (fonts, strips, margins, dots included) -- "
                        "e.g. 0.2 for a proportional 80%% -size miniature, not just a smaller canvas")
    p.add_argument("--no_dendrograms", action="store_true",
                   help="hide the dendrogram trees and reclaim their space for the heatmap; "
                        "clustering still determines row/column order, only the drawn tree goes away")
    p.add_argument("--label_pt_override", type=float, default=None,
                   help="hold tick-label font size fixed at this many points regardless of --scale "
                        "(geometry -- cell/strip/margin sizes -- still shrinks by scale as usual)")
    p.add_argument("--row_order_csv", type=Path, default=None,
                   help="single-column CSV/text file (no header) listing row (finch cluster) names "
                        "in the desired FIXED order, replacing hierarchical row clustering entirely -- "
                        "for matching this project's other, non-clustered heatmaps. Rows present in the "
                        "matrix but absent from this file are DROPPED, not appended (also serves as an "
                        "exclusion list, e.g. an edge-case cluster excluded elsewhere in the analysis)")
    main(**vars(p.parse_args()))
