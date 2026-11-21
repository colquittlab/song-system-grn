#!/usr/bin/env python
"""Zoom-in panel of MAFB-domain elements, laid out as the original submission's figure: zoom-ins side by side, small
accessibility tracks, each element in a grey band running through every track, the regulator motif hits in the element
listed below it in a regular grid (by position, filled column by column; not drawn at their coordinates), then
oscine and whole-tree phyloP.

Tracks are the lab hub's bigWigs, as in ucsc_dar_batch.py: accessibility (insertion, 10-bp bins, 20-bp smoothing) for
RA, C1H-1, PVALB-2, LGE-1 on a shared fixed axis (0-ACC_MAX, default 30) and per-base phyloP (oscine subtree, -3..3;
all 363, -4..6), each drawn as the mean per display bin.
Motifs: every regulator hit in the element from the curated-lobe scan (results/mafb_domain/element_motif_gain.tsv: the
12 SCENIC+ regulators, ARE-full, JASPAR AR) plus ARE-half (element_ar_sites.tsv, FIMO p < 1e-4); overlapping hits with
the same label are listed once. AR motifs (AREf, AREh, AR) are red; a hit with an oscine gain (oscine minus outgroup
presence >= 0.5, coverage floor met) is bold; hits in BOX are boxed.

    python zoom_panels.py                       # ZOOMS below -> results/ucsc/zoom_panels.{pdf,svg,png}
Overlapping hits of the AR matrices (ARE-full, ARE-half, JASPAR / SCENIC+ AR) are listed once.
"""
import math, os, re
import numpy as np, pandas as pd, pyBigWig
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                            "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results"
BASE = "https://colquitt-track-hub.pbsci.ucsc.edu/lonStrDom2/GCF_005870125.1"
W = f"{BASE}/scATAC/motor-pathway/wigs/insertion_bin10-smooth20"
ACC = [("Glut-RA", "#e31a1c", "RA"), ("Glut-Arco-1", "#fb9a99", "C1H-1"),
       ("GABA-4-1", "#4393c3", "PVALB-2"), ("GABA-1-1", "#92c5de", "LGE-1")]
PHY = [(f"{BASE}/comparative/subtree222_perbase_lsd2.bw", "#3f6b45", "#9a7b24", (-3, 3), "phyloP\noscine"),
       (f"{BASE}/comparative/wholetree_perbase_lsd2.bw", "#2166ac", "#b2182b", (-4, 6), "phyloP\navian")]
ACC_MAX = float(os.environ.get("ACC_MAX", "30"))
CHROM = "chr20"
# (label, [elements highlighted]); the window is the elements' span +/- FLANK
ZOOMS = [("1", ["RA:chr20:12714605-12715179"]),          # RA, oscine class; RA log2FC 3.5
         ("2", ["both:chr20:12732176-12733103"]),
         ("3", ["both:chr20:12945940-12947066#1", "both:chr20:12945940-12947066#2", "both:chr20:12945940-12947066#3"]),
         ("4", ["G4:chr20:13280687-13281096"])]          # PVALB-2, strong avian (+5.32), G4 log2FC 4.3
BOX = {12946573}                                        # motif hits (0-based start) boxed: the site shown in the alignment panel
FLANK = int(os.environ.get("FLANK", "400"))
LAB = {"full": "AREf", "half": "AREh", "jaspar": "AR", "AR": "AR"}
ARSET = {"AREf", "AREh", "AR"}
RED, INK, INK3, BAND = "#d7191c", "#171b22", "#6b7380", "#e6e6e6"
WIDTH_IN = float(os.environ.get("WIDTH_IN", "5.0"))
H_ACC, H_PHY, GAP = .15, .28, .03                # inches
FS, ROW_PT, COL_IN, MAXROWS, LEAD_PT = 5.0, 6.0, .27, 6, 3.0    # motif grid: font pt, row pitch pt, column width in, max rows per column
NBIN = 500


def motifs(elements):
    """-> {element: [(label, start, red, bold)]} sorted by position, same-label overlaps merged"""
    g = pd.read_csv(f"{R}/mafb_domain/element_motif_gain.tsv", sep="\t")
    h = pd.read_csv(f"{R}/mafb_domain/element_ar_sites.tsv", sep="\t")
    h = h[h.motif == "half"]
    out = {}
    for el in elements:
        rows = [(LAB.get(r.motif, r.motif), r.start, r.end, bool(r.floor_ok and r.gain >= .5))
                for r in g[g.element == el].itertuples()]
        rows += [("AREh", r.start, r.end, False) for r in h[h.element == el].itertuples()]
        rows.sort(key=lambda x: (x[1], x[0]))
        keep = []
        for lab, s, e, b in rows:                     # overlapping hits of one label, or of any AR matrices, listed once
            fam = lambda l: "AR" if l in ARSET else l
            dup = next((k for k in keep if fam(k[0]) == fam(lab) and s < k[2] and e > k[1]), None)
            if dup:
                # the label kept: the gained hit's, else ARE-full > AR > ARE-half
                pick = (dup[0] if dup[3] and not b else lab if b and not dup[3]
                        else min(dup[0], lab, key=["AREf", "AR", "AREh"].index) if fam(lab) == "AR" else lab)
                keep[keep.index(dup)] = (pick, min(dup[1], s), max(dup[2], e), dup[3] or b)
            else:
                keep.append((lab, s, e, b))
        out[el] = [(lab, s, lab in ARSET, b) for lab, s, e, b in keep]
    return out


def ext(el, T):
    r = T.loc[el]; return int(r.ext_start), int(r.ext_end)


def binned(url, a, b, n):
    v = pyBigWig.open(url).stats(CHROM, a, b, type="mean", nBins=n)
    return np.array([np.nan if x is None else x for x in v], float)


def windows():
    """-> [(label, start, end)] of the zoom-in windows (0-based, half-open), as drawn by main()"""
    T = pd.read_csv(f"{R}/mafb_domain/element_threshold_classes.tsv", sep="\t").set_index("element")
    out = []
    for lab, els in ZOOMS:
        sp = [ext(e, T) for e in els]
        out.append((lab, min(s for s, _ in sp) - FLANK, max(e for _, e in sp) + FLANK))
    return out


def main():
    T = pd.read_csv(f"{R}/mafb_domain/element_threshold_classes.tsv", sep="\t").set_index("element")
    M = motifs([e for _, els in ZOOMS for e in els])
    shape = lambda n: (math.ceil(n / MAXROWS), math.ceil(n / math.ceil(n / MAXROWS)))   # (columns, rows), rows balanced
    nrows = max(shape(len(v))[1] for v in M.values() if v)
    h_grid = (nrows * ROW_PT + LEAD_PT + 4) / 72
    tracks = [("acc", i) for i in range(len(ACC))] + [("grid", 0)] + [("phy", i) for i in range(len(PHY))]
    hts = [H_ACC] * len(ACC) + [h_grid] + [H_PHY] * len(PHY)
    top_in, bot_in, left_in, right_in, gap_x = .3, .22, .52, .03, .1
    H = top_in + sum(hts) + GAP * (len(hts) - 1) + bot_in
    fig = plt.figure(figsize=(WIDTH_IN, H))
    pw = (WIDTH_IN - left_in - right_in - gap_x * (len(ZOOMS) - 1)) / len(ZOOMS)
    for z, (lab, els) in enumerate(ZOOMS):
        spans = [ext(e, T) for e in els]
        a = min(s for s, _ in spans) - FLANK; b = max(e for _, e in spans) + FLANK
        x0 = left_in + z * (pw + gap_x)
        acc = [binned(f"{W}/{f}-insertion-bin10-smooth20.bw", a, b, NBIN) for f, _, _ in ACC]
        phy = [binned(u, a, b, NBIN) for u, *_ in PHY]
        xs = np.linspace(a, b, NBIN + 1)
        y = H - top_in
        axes = []
        for (kind, i), ht in zip(tracks, hts):
            y -= ht
            ax = fig.add_axes([x0 / WIDTH_IN, y / H, pw / WIDTH_IN, ht / H]); axes.append((kind, i, ax))
            y -= GAP
            ax.set_xlim(a, b); ax.axis("off")
            for s, e in spans:
                ax.add_patch(Rectangle((s, -1e3), e - s, 2e3, color=BAND, lw=0, zorder=0))
            if kind == "acc":
                v = np.clip(np.nan_to_num(acc[i]), 0, ACC_MAX)
                ax.stairs(v, xs, fill=True, color=ACC[i][1], lw=0, zorder=2)
                ax.set_ylim(0, ACC_MAX)
                ax.axhline(0, color=ACC[i][1], lw=.3, zorder=2)
            elif kind == "phy":
                url, cp, cn, (lo, hi), name = PHY[i]
                v = np.clip(np.nan_to_num(phy[i]), lo, hi)
                ax.stairs(np.maximum(v, 0), xs, fill=True, color=cp, lw=0, zorder=2)
                ax.stairs(np.minimum(v, 0), xs, fill=True, color=cn, lw=0, zorder=2)
                ax.set_ylim(lo, hi); ax.axhline(0, color="#9aa0a6", lw=.3, zorder=1)
            else:                                     # motif grid under each element, placed left to right without overlap
                ax.set_ylim(0, 1)
                bx = ax.get_position(); px_per_bp = pw / (b - a)
                placed_right = -1e9
                for (s, e), el in zip(spans, els):
                    ms = M[el]
                    if not ms: continue
                    ncol, nr = shape(len(ms))
                    wbp = ncol * COL_IN / px_per_bp
                    left = max((s + e) / 2 - wbp / 2, placed_right, a)
                    left = min(left, b - wbp)
                    placed_right = left + wbp
                    lead = LEAD_PT / 72 / h_grid             # room above the rows for a leader line
                    gc = left + wbp / 2
                    if abs(gc - (s + e) / 2) > max((e - s) / 2, .25 * COL_IN / px_per_bp):   # grid moved off its band
                        ax.plot([(s + e) / 2, (s + e) / 2, gc], [1, 1 - lead * .45, 1 - lead * .45], color=INK3, lw=.35,
                                zorder=3, solid_joinstyle="miter")
                    for k, (name, pos, red, bold) in enumerate(ms):
                        c, r = divmod(k, nr)
                        ax.text(left + (c + .5) * COL_IN / px_per_bp, 1 - lead - (r + .5) * ROW_PT / 72 / h_grid, name,
                                ha="center", va="center", fontsize=FS, color=RED if red else INK,
                                fontweight="bold" if bold else "normal", zorder=3,
                                bbox=dict(boxstyle="square,pad=0.18", fc="none", ec=INK, lw=.5) if pos in BOX else None)
        # header: zoom label and a scale bar
        ax0 = axes[0][2]
        sb = 1000 if b - a >= 3000 else 500 if b - a >= 1600 else 250
        ax0.plot([a + (b - a) * .02, a + (b - a) * .02 + sb], [ACC_MAX * 1.55] * 2, color=INK, lw=.6, clip_on=False)
        ax0.text(a + (b - a) * .02 + sb / 2, ACC_MAX * 1.62, f"{sb // 1000} kb" if sb >= 1000 else f"{sb} bp",
                 ha="center", va="bottom", fontsize=FS, color=INK, clip_on=False)
        ax0.text((a + b) / 2, ACC_MAX * 1.62, lab, ha="center", va="bottom", fontsize=7, fontweight="bold", color=INK,
                 clip_on=False)
        print(f"  zoom-in {lab}: {CHROM}:{a + 1:,}-{b:,}")
        if z == 0:                                    # track names and y limits, left of the first zoom-in
            for kind, i, ax in axes:
                if kind == "grid": continue
                if kind == "acc": name, (lo, hi) = ACC[i][2], (0, ACC_MAX)
                else: name, (lo, hi) = PHY[i][4], PHY[i][3]
                ax.text(-.11, .5, name, transform=ax.transAxes, ha="right", va="center", fontsize=FS + .5, color=INK,
                        linespacing=.95)
                for v_, va in ((hi, "top"), (lo, "bottom")):
                    ax.text(-.015, (v_ - lo) / (hi - lo), f"{v_:g}", transform=ax.transAxes, ha="right", va=va,
                            fontsize=FS - .8, color=INK3)
    fig.text(left_in / WIDTH_IN, .03 / H, "motif hits in each element, listed by position (not to scale); overlapping AR "
             "hits listed once\nred: AR motifs (AREf ARE-full, AREh ARE-half, AR JASPAR / SCENIC+);  bold: site gained "
             "in oscines (oscine minus outgroup presence ≥ 0.5)", fontsize=FS - .5, color=INK3, va="bottom", linespacing=1.2)
    os.makedirs(f"{R}/ucsc", exist_ok=True)
    stem = f"{R}/ucsc/zoom_panels"
    for ext_ in ("pdf", "svg", "png"):
        fig.savefig(f"{stem}.{ext_}", dpi=400)
    print(f"  wrote {stem}.pdf / .svg / .png  ({WIDTH_IN:g} x {H:.2f} in)")
    for el, v in M.items():
        print(f"    {el}: " + ", ".join(f"{n}{'*' if b else ''}" for n, _, _, b in v))


if __name__ == "__main__":
    main()
