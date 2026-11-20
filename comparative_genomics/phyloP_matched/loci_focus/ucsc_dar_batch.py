#!/usr/bin/env python
"""UCSC Genome Browser PDFs of each MAFB-domain DAR peak, from the lab hub's own bigWigs.

For every SnapATAC2 peak in results/mafb_domain/element_summaries.tsv (RA, G4 and shared elements;
one view per peak, lobes not split), the browser draws a window of constant width (default 10 kb)
centred on the peak's midpoint, with only these tracks, in
this order, as custom tracks pointing at the hub's bigWig files (so order, colours and scales are set
here rather than by the hub's trackDb):
  accessibility (insertion, 10-bp bins, 20-bp smoothing), labelled RA (Glut-RA), C1H-1 (Glut-Arco-1),
  PVALB-2 (GABA-4-1) and LGE-1 (GABA-1-1)
  oscine-lineage phyloP (subtree birdAnc222) and whole-tree phyloP (363 birds), per base
All four accessibility tracks share a fixed y axis, 0 to ACC_MAX (default 30; values above are
clipped). With ACC_MAX=pair, Glut-RA/Glut-Arco-1 and GABA-4-1/GABA-1-1 instead each share 0 to 1.05 x the
pair's maximum over the view. phyloP limits follow the hub (oscine -3:3, whole tree -4:6). Vertical guidelines and the
descriptions above tracks are off; tracks are HEIGHT px tall (45, 1.5 x the hub default); the peak is highlighted in grey. The UCSC API key is read from
~/.ucsc_api_key (never stored here).

    ACC_MAX=30 python ucsc_dar_batch.py [window_bp=10000] [first N peaks]
    python ucsc_dar_batch.py locus SIX2 [flank_bp=10000]   # one gene +/- flank, peak-batch tracks + gene models
    python ucsc_dar_batch.py region NAME chr:start-end [chr:start-end ...]   # one element (grey) + marked sites (darker),
                                          # 10-kb and 2-kb windows, peak-batch tracks
    python ucsc_dar_batch.py ar_sites [win_bp=2000]   # every significant oscine AR-site gain in the curated lobes,
                                          # one view each with y-axis limits -> results/ucsc/ar_sites_win2k_acc30/
    python ucsc_dar_batch.py overview     # the whole domain (chr20:12,551,036-13,547,533): accessibility (autoscaled,
                                          # 30 px) and gene models, no highlights -> results/ucsc/mafb_domain_overview.pdf
Writes results/ucsc/dar_peaks_win<W>k_acc<max>/<nn>_<class>_chr20_<start>-<end>.pdf and the merged
dar_peaks_win<W>k_acc<max>.pdf (all pages, in genome order).
"""
import os, re, subprocess, sys, time, urllib.parse, urllib.request, http.cookiejar
import numpy as np, pandas as pd, pyBigWig

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results"
BASE = "https://colquitt-track-hub.pbsci.ucsc.edu/lonStrDom2"
HUB, GENOME, UCSC = f"{BASE}/hub.txt", "GCF_005870125.1", "https://genome.ucsc.edu"
W = f"{BASE}/{GENOME}/scATAC/motor-pathway/wigs/insertion_bin10-smooth20"
# (cluster file, colour, scale pair, label): RA red / C1H-1 light red as on the hub; PVALB-2 (GABA-4-1) blue
# #4393C3 as in the phyloP figures, LGE-1 (GABA-1-1) the lighter #92C5DE
ACC = [("Glut-RA", "227,26,28", "RA", "RA"), ("Glut-Arco-1", "251,154,153", "RA", "C1H-1"),
       ("GABA-4-1", "67,147,195", "G4", "PVALB-2"), ("GABA-1-1", "146,197,222", "G4", "LGE-1")]
PHY = [("phyloP oscine", f"{BASE}/{GENOME}/comparative/subtree222_perbase_lsd2.bw", "63,107,69", "154,123,36", "-3:3"),
       ("phyloP all 363", f"{BASE}/{GENOME}/comparative/wholetree_perbase_lsd2.bw", "33,102,172", "178,24,43", "-4:6")]
PIX = 1200
HEIGHT = 45
OVERVIEW_PIX = int(os.environ.get("OVERVIEW_PIX", "1600"))  # image width for the domain overview: a full browser width,
OVERVIEW_TEXT = os.environ.get("OVERVIEW_TEXT", "18")         # so 2-px smoothing spans what it does on screen; larger browser
                                                              # text keeps labels readable once scaled to OVERVIEW_IN
OVERVIEW_HEIGHT = int(os.environ.get("OVERVIEW_HEIGHT", "44")) # accessibility track height (px) at that width
OVERVIEW_END = 13545500                                       # track limits come from the domain up to ~2 kb short of the TOP1
                                                              # TSS (13,547,533): its promoter peak (RA 93, PVALB-2 174) is clipped
OVERVIEW_VIEW_END = 13549500                                  # the view reaches just past the TOP1 TSS, so its gene model shows                                          # track height (px): 1.5 x the hub default of 30


GENES = f"{BASE}/{GENOME}/genes/GCF_005870125.1_lonStrDom2_genomic_stringtie_FPKMthresh0.5_FPKMmaxfrac0.1_minhits5_slim.bb"


def gene_track(chrom, a, b):
    """inline BED12 of the hub's gene models over [a, b): the longest transcript per gene, named by symbol (name2)"""
    bb = pyBigWig.open(GENES); rows = {}
    for s, e, rest in bb.entries(native(chrom), a, b) or []:
        f = rest.split("\t"); sym = f[9] or f[0]
        if sym not in rows or e - s > rows[sym][1] - rows[sym][0]: rows[sym] = (s, e, f)
    out = ['track name="genes" description="NCBI gene models (3\' UTR extended), longest transcript per gene" '
           'visibility=pack color=60,60,60']
    for sym, (s, e, f) in sorted(rows.items(), key=lambda kv: kv[1][0]):
        out.append("\t".join(map(str, [chrom, s, e, sym, 0, f[2], f[3], f[4], 0, f[6], f[7], f[8]])))
    return "\n".join(out)


NATIVE = {"chr20": "NC_042585.1"}


def native(chrom):
    """chrN -> the RefSeq name the gene bigBed uses, matched by sequence length (cached in NATIVE)"""
    if chrom not in NATIVE:
        L = pyBigWig.open(PHY[1][1]).chroms()[chrom]
        NATIVE[chrom] = next(c for c, l in pyBigWig.open(GENES).chroms().items() if l == L)
    return NATIVE[chrom]


def find_gene(sym):
    """-> (chrN, start, end, strand) of the longest gene model named sym"""
    bb = pyBigWig.open(GENES); rev = {}; best = None
    bw = pyBigWig.open(PHY[1][1]).chroms(); byl = {l: c for c, l in bw.items()}
    for c, L in bb.chroms().items():
        for s, e, rest in bb.entries(c, 0, L) or []:
            f = rest.split("\t")
            if (f[9] or f[0]) == sym and (best is None or e - s > best[2] - best[1]):
                best = (byl[L], s, e, f[2]); NATIVE[byl[L]] = c
    if best is None: raise SystemExit(f"no gene model named {sym}")
    return best


def locus(key, sym, flank=10000):
    """one gene with `flank` bp either side, in the DAR-peak style (four accessibility tracks on a shared fixed axis,
    oscine and whole-tree phyloP, HEIGHT px) plus the hub's gene models, no highlight. Renders at ACC_MAX (default 30,
    as the peak batch) and, when the view's accessibility exceeds it, again at the view's maximum rounded up to 10.
    Writes results/ucsc/<sym>_flank<k>k_acc<max>_browser.pdf (the browser's PDF) and .svg/.pdf (editable text, labels
    black and right-aligned, cropped; LOCUS_IN wide, default 6.8 in)."""
    chrom, s, e, strand = find_gene(sym); a, b = max(0, s - flank), e + flank
    top = max(pyBigWig.open(f"{W}/{n}-insertion-bin10-smooth20.bw").stats(chrom, a, b, type="max", exact=True)[0] or 0 for n, *_ in ACC)
    print(f"  {sym}: {chrom}:{s + 1:,}-{e:,} ({strand}); view {chrom}:{a + 1:,}-{b:,}; accessibility max {top:.1f}")
    scales = [float(os.environ.get("ACC_MAX", "30"))]
    if top > scales[0]: scales.append(float(np.ceil(top / 10) * 10))
    for sc in scales:
        stem = f"{R}/ucsc/{sym}_flank{flank / 1000:g}k_acc{sc:g}"
        render(key, chrom, a, b, None, None, {"RA": sc, "G4": sc}, f"{stem}_browser.pdf", highlight="",
               extra=gene_track(chrom, a, b))
        w, h = to_figure(f"{stem}_browser.pdf", stem, float(os.environ.get("LOCUS_IN", "6.8")))
        print(f"    0-{sc:g}: wrote {stem}_browser.pdf and {stem}.svg/.pdf ({w:g} x {h:.2f} in)")
P2G = f"{BASE}/{GENOME}/scATAC/motor-pathway/peaks2genes/hybrid_maxdist1Mb/peak2geneLinks_ra_corCutOff-{{cut}}.bb"


def p2g_track(chrom, a, b, gene="MAFB", cut="0.45"):
    """inline interact track of the hub's 1-Mb RA + C1H-1 peak-to-gene links (r >= cut) to `gene`, peak in [a, b)"""
    bb = pyBigWig.open(P2G.format(cut=cut)); rows = []
    for s, e, rest in bb.entries(chrom, a, b):
        f = rest.split("\t")                      # name score value exp color src* tgt*
        g, peak = f[0].split("<>")
        ps = int(peak.split("_")[1].split("-")[0])
        if g == gene and a <= ps < b:            # interact custom tracks need the native sequence name, not the chr20 alias
            t = min(1.0, max(0.0, (float(f[2]) - float(cut)) / (1 - float(cut))))   # r: cut -> light grey, 1 -> near black
            lo, hi = (210, 210, 210), (25, 25, 25)
            f[4] = "#%02x%02x%02x" % tuple(round(l + t * (h - l)) for l, h in zip(lo, hi))
            f[0] = f[8] = f[13] = "."                 # no labels on the arcs
            rows.append([NATIVE.get(chrom, chrom), s, e] + [NATIVE.get(x, x) for x in f])
    out = [f'track type=interact name="RA + C1H-1 P2G r >= {cut}" description="RA + C1H-1 peak-to-gene links to {gene}, '
           f'r >= {cut}, max distance 1 Mb; grey, light (r = {cut}) to black (r = 1)" visibility=full maxHeightPixels=200:60:20 interactUp=true']
    out += ["\t".join(map(str, r)) for r in rows]
    return "\n".join(out), len(rows)


def custom_text(lim, phy=True, height=None, auto=False, extra="", pre="", smooth=None, tlim=None):
    height = height or HEIGHT
    # per-pixel mean, then smoothed: the default here is mean+whiskers, on which smoothingWindow has almost no effect
    sm = ("windowingFunction=mean " if (smooth or tlim) else "") + (f"smoothingWindow={smooth} " if smooth else "")
    t = []
    for name, col, pair, lab in ACC:
        t.append(f'track type=bigWig name="{lab}" bigDataUrl={W}/{name}-insertion-bin10-smooth20.bw color={col} '
                 + (f"autoScale=off viewLimits={tlim[name][0]:.3g}:{tlim[name][1]:.3g} alwaysZero=off " if tlim
                    else "autoScale=on " if auto else f"autoScale=off viewLimits=0:{lim[pair]:g} ") + sm
                 + ("" if tlim else "alwaysZero=on ") + f'maxHeightPixels=128:{height}:11 visibility=full')
    for name, url, col, alt, vl in (PHY if phy else []):
        t.append(f'track type=bigWig name="{name}" bigDataUrl={url} color={col} altColor={alt} autoScale=off '
                 f'viewLimits={vl} alwaysZero=on windowingFunction=mean maxHeightPixels=128:{HEIGHT}:11 visibility=full')
    return (pre + "\n" if pre else "") + "\n".join(t) + ("\n" + extra if extra else "")


def render(key, chrom, a, b, hs, he, lim, out, phy=True, highlight=None, height=None, auto=False, show=None, extra="", pre="",
           pix=None, smooth=None, text=None, tlim=None, axis=False):
    """axis: centre labels on, which is the only way the browser draws each track's y-axis limits (the centre text
    itself is removed later by tidy(axis=True)).
    show: hub tracks to keep visible as {track: visibility}; highlight "" = none"""
    """two requests in one session: load the custom tracks (their ids are assigned by the server), then
    hide every track and show only those ids -- hideTracks would otherwise hide the custom tracks too"""
    cj = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    base = {"hubUrl": HUB, "genome": GENOME, "position": f"{chrom}:{a + 1}-{b}", "pix": pix or PIX, "apiKey": key,
            **({"textSize": text} if text else {})}
    post = lambda q: cj.open(f"{UCSC}/cgi-bin/hgTracks", urllib.parse.urlencode(q).encode(), timeout=300).read().decode("utf-8", "replace")
    for i in range(3):
        try:
            page = post({**base, "hgt.customText": custom_text(lim, phy, height, auto, extra, pre, smooth, tlim)})
            ids = list(dict.fromkeys(re.findall(r"ct_[A-Za-z0-9]+_[0-9]+", page)))
            if len(ids) < len(ACC) + (len(PHY) if phy else 0) + bool(extra) + bool(pre): raise RuntimeError(f"custom tracks not loaded ({ids})")
            q = {**base, "hideTracks": 1, **{t: ("pack" if t.startswith("ct_genes_") else "full") for t in ids}, **(show or {}),
                 "guidelines": "off", **({} if axis else {"centerLabels": "off"}),
                 "hgt.psOutput": "on"}
            if highlight != "": q["highlight"] = highlight or f"{chrom}:{hs + 1}-{he}#e6e6e6"
            page = post(q)
            p = re.search(r"\.\./trash/hgt/hgt_genome_[^'\"]*\.pdf", page)
            if not p: raise RuntimeError("no PDF link in response")
            open(out, "wb").write(cj.open(f"{UCSC}/{p.group(0)[3:]}", timeout=120).read()); return
        except Exception:
            if i == 2: raise
            time.sleep(10 * (i + 1))


CLEAN = "/private/home/bcolquit/group/song-system-grn/scripts/browser_pdf_to_svg.py"
PY_LXML = "/private/home/bcolquit/miniforge3/envs/crested/bin/python"     # an interpreter with lxml for CLEAN


def tidy(t, pad=3.0, axis=False, numw=9.0):
    """drop the browser's pale-red border at the left of the data area, the end markers of interact arcs and the magenta
    clipping marks; set the
    left-column track labels in black, right-aligned flush (pad units) against the data area"""
    import re
    m = re.search(r'<path\b[^>]*fill:#ffb4b4[^>]*>', t)
    x0 = 139.5
    if m:                                             # the border: a full-height 6.75-unit bar; data starts at its right edge
        d = re.search(r'\sd="m ([0-9.]+),[0-9.]+ h ([0-9.]+)', m.group(0)); sc = re.search(r'matrix\(([0-9.]+)', m.group(0))
        if d and sc: x0 = (float(d.group(1)) + float(d.group(2))) * float(sc.group(1))
        t = t.replace(m.group(0), "", 1)
    # interact end markers: small filled boxes in a thin band above the arcs' baseline (raw, y-up coordinates)
    def ybox(d):
        x = y = 0.0; ys = []
        for c, a in re.findall(r'([MmLlHhVvZz])([^MmLlHhVvZzCcSsQqTtAa]*)', d):
            v = [float(q) for q in re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', a)]
            if c in "Mm" or c in "Ll":
                for i in range(0, len(v) - 1, 2):
                    x, y = (v[i], v[i + 1]) if c.isupper() else (x + v[i], y + v[i + 1]); ys.append(y)
            elif c in "Hh": x = v[-1] if c == "H" else x + sum(v)
            elif c in "Vv":
                for q in v: y = q if c == "V" else y + q; ys.append(y)
        return (min(ys), max(ys)) if ys else None
    base = [float(re.search(r'\sd="m [0-9.]+,([0-9.]+) c', el).group(1)) for el in re.findall(r'<path\b[^>]*/>', t)
            if re.search(r'\sd="m [0-9.]+,[0-9.]+ c', el) and re.search(r'stroke:#[0-9a-f]{6}', el)]   # arc starts
    if base:
        b0 = min(base)
        def drop(mm):
            el = mm.group(0)
            if "stroke:none" not in el or re.search(r'\sd="[^"]*[CcSsQqAa]', el): return el
            d = re.search(r'\sd="([^"]*)"', el); yb = ybox(d.group(1)) if d else None
            return "" if yb and yb[0] >= b0 - 1 and yb[1] <= b0 + 60 else el
        t = re.sub(r'<path\b[^>]*/>', drop, t)
    # the browser's magenta marks for values clipped by the view limits (below a raised floor)
    t = re.sub(r'<path\b[^>]*fill:#ff00ff[^>]*/>', "", t)
    if axis:
        return axis_labels(t, x0, pad, numw)
    def fix(mm):                                      # texts in the left label column (translate x < 20)
        el = mm.group(0); tr = re.search(r'transform="matrix\(([0-9.]+),0,0,([0-9.]+),([-0-9.]+),([-0-9.]+)\)"', el)
        if not tr or float(tr.group(3)) >= 20: return el
        el = el.replace(tr.group(0), f'transform="matrix({tr.group(1)},0,0,{tr.group(2)},{x0 - pad:.3f},{tr.group(4)})"')
        return re.sub(r'fill:#[0-9a-fA-F]{6}', 'fill:#000000', el).replace('writing-mode:lr-tb;', 'writing-mode:lr-tb;text-anchor:end;')
    return re.sub(r'<text\b.*?</text>', fix, t, flags=re.S)


def axis_labels(t, x0, pad, numw):
    """labels when the browser drew y-axis limits (centre labels on): drop the centre labels (track names over the data
    area); set each track's limits right-aligned flush against the data area, without the browser's "_" ticks, and its
    name right-aligned just left of them, centred between them; everything black"""
    import re
    els = re.findall(r'<text\b.*?</text>', t, re.S)
    tsp = lambda el: re.findall(r'<tspan\b([^>]*)>(.*?)</tspan>', el, re.S)
    TR = r'transform="matrix\(([-0-9.]+),0,0,([-0-9.]+),([-0-9.]+),([-0-9.]+)\)"'
    names = set()
    for el in els:                                    # the left label column: a name, then limits ending in "_"
        tr = re.search(TR, el)
        if tr and float(tr.group(3)) < 20:
            names |= {c.strip() for _, c in tsp(el) if not c.strip().endswith("_")}
    def one(attr, text, x, y, anchor):
        st = re.search(r'style="([^"]*)"', attr).group(1)
        st = re.sub(r'fill:#[0-9a-fA-F]{6}', 'fill:#000000', st).replace('writing-mode:lr-tb;', f'writing-mode:lr-tb;text-anchor:{anchor};')
        return f'<tspan style="{st}" x="{x:.3f}" y="{y:.3f}">{text}</tspan>'
    def fix(mm):
        el = mm.group(0); tr = re.search(TR, el)
        if not tr: return el
        sx, sy, tx, ty = (float(v) for v in tr.groups())
        parts = tsp(el)
        if tx >= x0:                                  # data area: the centre labels go
            return "" if "".join(c for _, c in parts).strip() in names else el
        lim = [(a, c) for a, c in parts if c.strip().endswith("_")]
        nm = [(a, c) for a, c in parts if not c.strip().endswith("_")]
        if tx >= 20 and nm: return el                 # Scale / chromosome / ruler text
        out = []
        for a, c in lim:
            y = float(re.search(r'\sy="([-0-9.]+)"', a).group(1))
            # the browser's "_" marks the edge: the top limit hangs below it (cap height ~0.72 em), the bottom sits on it
            out.append(one(a, c.strip().rstrip("_").strip(), 0, y + (4.3 if y < 0 else 0.0), "end"))
        for a, c in nm:                               # name centred on the track (top limit -8.28, bottom +11.97)
            out.append(one(a, c.strip(), -numw, (-8.28 + 11.97) / 2 + 2.0, "end"))
        return (f'<text xml:space="preserve" transform="matrix({sx},0,0,{sy},{x0 - pad:.3f},{ty})">'
                + "".join(out) + "</text>")
    return re.sub(r'<text\b.*?</text>', fix, t, flags=re.S)


def to_figure(pdf, stem, width_in, axis=False):
    """browser PDF -> SVG with editable text (Inkscape, fonts kept) -> merged paths (CLEAN) -> tidied -> cropped to the
    drawing -> scaled to width_in
    inches -> PDF (Inkscape). Writes <stem>.svg and <stem>.pdf."""
    import re, tempfile
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        sh = lambda c: subprocess.run(["bash", "-lc", c], check=True, capture_output=True)
        sh(f"inkscape '{pdf}' --pdf-font-strategy=keep --export-type=svg --export-plain-svg --export-filename='{d}/raw.svg'")
        subprocess.run([PY_LXML, CLEAN, f"{d}/raw.svg", "-o", f"{d}/clean.svg"], check=True, capture_output=True)
        open(f"{d}/tidy.svg", "w").write(tidy(open(f"{d}/clean.svg").read(), axis=axis))
        # crop to the drawn content (the label column is wider than the right-aligned labels need)
        sh(f"inkscape '{d}/tidy.svg' --export-area-drawing --export-plain-svg --export-type=svg --export-filename='{d}/crop.svg'")
        t = open(f"{d}/crop.svg").read()
        m = re.search(r'<svg[^>]*?\swidth="([0-9.]+)"[^>]*?\sheight="([0-9.]+)"', t, re.S)
        w, h = float(m.group(1)), float(m.group(2))
        t = t.replace(m.group(0), m.group(0).replace(f'width="{m.group(1)}"', f'width="{width_in}in"')
                                             .replace(f'height="{m.group(2)}"', f'height="{width_in * h / w:.4f}in"'), 1)
        open(f"{stem}.svg", "w").write(t)
        sh(f"inkscape '{stem}.svg' --export-type=pdf --export-filename='{stem}.pdf'")
    return width_in, width_in * h / w


def track_limits(chrom, a, b, pix, smooth, floor_q=None, lim_end=None):
    """per-track view limits mimicking autoscale with a raised floor: top = the highest per-pixel mean (smoothed as
    drawn), bottom = that track's floor_q-th percentile pixel (OVERVIEW_FLOOR, default the median)"""
    floor_q = float(os.environ.get("OVERVIEW_FLOOR", "50")) if floor_q is None else floor_q
    out = {}
    for name, *_ in ACC:
        v = np.nan_to_num(np.array(pyBigWig.open(f"{W}/{name}-insertion-bin10-smooth20.bw").stats(chrom, a, b, type="mean", nBins=pix), float))
        if smooth: v = np.convolve(v, np.ones(smooth) / smooth, "same")
        if lim_end: v = v[: int(len(v) * (lim_end - a) / (b - a))]   # limits from the pixels before lim_end only
        out[name] = (float(np.percentile(v, floor_q)), float(v.max()))
    return out


def overview(key, out, smooth=2):
    """the whole DHX35-TOP1 domain: RA + C1H-1 peak-to-gene links to MAFB (1-Mb run, r >= 0.45) above the
    accessibility tracks (per-pixel mean, smoothed over `smooth` pixels or not at all; each track scaled from its
    OVERVIEW_FLOOR-th percentile pixel (default 50, the median) to its maximum, lifting the background noise), then the hub's gene models. Track labels
    are set black and right-aligned against the tracks, and the browser's red left border and the arcs' end markers are
    removed (tidy). The browser PDF (always 7.5 in wide) is kept as <out>_browser.pdf; the figure is OVERVIEW_IN wide."""
    a, b = 12551035, OVERVIEW_VIEW_END               # domain start (DHX35) to just past the TOP1 TSS
    pre, n = p2g_track("chr20", a, OVERVIEW_END)
    tlim = track_limits("chr20", a, b, OVERVIEW_PIX, smooth, lim_end=OVERVIEW_END)   # the TOP1 promoter peak is clipped
    raw = out.replace(".pdf", "_browser.pdf")
    render(key, "chr20", a, b, None, None, {}, raw, phy=False, highlight="", height=OVERVIEW_HEIGHT, auto=True,
           extra=gene_track("chr20", a, b), pre=pre, pix=OVERVIEW_PIX, smooth=smooth, text=OVERVIEW_TEXT, tlim=tlim)
    w, h = to_figure(raw, out[:-4], float(os.environ.get("OVERVIEW_IN", "6.8")))
    print(f"  {'smoothed ' + str(smooth) + ' px' if smooth else 'unsmoothed'}: limits " + ", ".join(f"{k} {lo:.2f}-{hi:.1f}" for k, (lo, hi) in tlim.items())
          + f"; {n} links; wrote {out[:-4]}.svg/.pdf ({w:g} x {h:.2f} in)")


def region(key, name, chrom, s, e, marks=(), wins=(10000, 2000)):
    """one element [s, e) in the DAR-peak style (fixed 0-ACC_MAX accessibility, phyloP, HEIGHT px), windows of constant
    width centred on it; the element highlighted in grey, each (start, end) in marks in a darker grey. Writes
    results/ucsc/<name>_win<k>k_acc<max>_browser.pdf and the tidied .svg/.pdf (REGION_IN wide, default 6.8 in)."""
    sc = float(os.environ.get("ACC_MAX", "30"))
    hl = "|".join([f"{chrom}:{s + 1}-{e}#e6e6e6"] + [f"{chrom}:{a + 1}-{b}#bdbdbd" for a, b in marks])
    for win in wins:
        mid = (s + e) // 2; a = max(0, mid - win // 2); b = a + win
        stem = f"{R}/ucsc/{name}_win{win / 1000:g}k_acc{sc:g}"
        render(key, chrom, a, b, s, e, {"RA": sc, "G4": sc}, f"{stem}_browser.pdf", highlight=hl)
        w, h = to_figure(f"{stem}_browser.pdf", stem, float(os.environ.get("REGION_IN", "6.8")))
        print(f"  {name}: view {chrom}:{a + 1:,}-{b:,}; wrote {stem}_browser.pdf and {stem}.svg/.pdf ({w:g} x {h:.2f} in)")


def ar_sites(key, win=2000, q_max=0.05):
    """every AR site (ARE-full, JASPAR AR, SCENIC+ AR) in the curated MAFB-domain lobes whose oscine presence exceeds
    the outgroup's (one-sided Fisher, BH q < q_max over every site clearing the coverage floor, all matrices;
    results/mafb_domain/element_motif_gain.tsv). Overlapping hits of different AR matrices are one locus and one view:
    win bp centred on it, its lobe in grey, the site in darker grey, y-axis limits shown (axis=True). Writes
    results/ucsc/ar_sites_win<k>k_acc<max>/<nn>_<class>_chr20_<start>-<end>_<matrices>{_browser.pdf,.svg,.pdf}, the
    merged PDF of the figures and index.tsv (one row per matrix hit)."""
    d = pd.read_csv(f"{R}/mafb_domain/element_motif_gain.tsv", sep="\t")
    d = d[d.floor_ok].copy()
    p = d.fisher_p.values; o = np.argsort(p); q = np.empty(len(p))
    q[o] = np.minimum(np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1], 1)
    d["q"] = q
    a = d[d.motif.isin(["full", "jaspar", "AR"]) & (d.q < q_max)].sort_values(["start", "end"])
    loci = []                                         # overlapping hits merged
    for r in a.itertuples(index=False):
        if loci and r.start < loci[-1]["end"] and r.element == loci[-1]["element"]:
            L = loci[-1]; L["end"] = max(L["end"], r.end); L["rows"].append(r)
        else:
            loci.append(dict(element=r.element, cls=r.cls, chrom=r.chrom, start=r.start, end=r.end, rows=[r]))
    el = pd.read_csv(f"{R}/mafb_domain/element_threshold_classes.tsv", sep="\t").set_index("element")
    sc = float(os.environ.get("ACC_MAX", "30"))
    out = f"{R}/ucsc/ar_sites_win{win / 1000:g}k_acc{sc:g}"; os.makedirs(out, exist_ok=True)
    MN = {"full": "AREfull", "jaspar": "JASPAR", "AR": "SCENIC-AR"}
    idx, figs = [], []
    for i, L in enumerate(loci, 1):
        e = el.loc[L["element"]]; cls = {"both": "shared"}.get(L["cls"], L["cls"])
        mats = "+".join(MN[m] for m in dict.fromkeys(r.motif for r in L["rows"]))
        stem = f"{out}/{i:02d}_{cls}_{L['chrom']}_{L['start'] + 1}-{L['end']}_{mats}"
        mid = (L["start"] + L["end"]) // 2; a0 = max(0, mid - win // 2); b0 = a0 + win
        hl = f"{L['chrom']}:{int(e.ext_start) + 1}-{int(e.ext_end)}#e6e6e6|{L['chrom']}:{L['start'] + 1}-{L['end']}#bdbdbd"
        render(key, L["chrom"], a0, b0, None, None, {"RA": sc, "G4": sc}, f"{stem}_browser.pdf", highlight=hl, axis=True)
        to_figure(f"{stem}_browser.pdf", stem, float(os.environ.get("REGION_IN", "6.8")), axis=True); figs.append(f"{stem}.pdf")
        for r in L["rows"]:
            idx.append(dict(n=i, file=os.path.basename(stem) + ".pdf", element=L["element"], dar_class=cls,
                            phylop_class=e.thr_class, depth_ok=r.depth_ok, matrix=r.motif, site=f"{r.chrom}:{r.start + 1}-{r.end}",
                            osc=f"{r.p_osc}/{r.n_osc}", out=f"{r.p_out}/{r.n_out}", gain=round(r.gain, 3),
                            fisher_p=r.fisher_p, q=r.q, tier=r.tier, principal=r.principal))
        print(f"  {i:2d}/{len(loci)} {cls:6s} {L['chrom']}:{L['start'] + 1:,}-{L['end']:,} {mats:24s} "
              + ", ".join(f"{r.motif} {r.gain:+.2f} q={r.q:.1e} {r.tier}" for r in L["rows"]), flush=True)
        time.sleep(2)
    pd.DataFrame(idx).to_csv(f"{out}/index.tsv", sep="\t", index=False, float_format="%.3g")
    subprocess.run(["pdfunite", *figs, f"{out}.pdf"], check=True)
    print(f"  {len(a)} AR hits at q < {q_max} in {len(loci)} loci; wrote {out}/ (index.tsv) and {out}.pdf")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "ar_sites":
        key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
        ar_sites(key, int(sys.argv[2]) if len(sys.argv) > 2 else 2000); return
    if len(sys.argv) > 1 and sys.argv[1] == "region":
        key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
        os.makedirs(f"{R}/ucsc", exist_ok=True)
        c, rng = sys.argv[3].split(":"); s, e = (int(x.replace(",", "")) for x in rng.split("-"))
        marks = [tuple(int(x.replace(",", "")) for x in m.split(":")[1].split("-")) for m in sys.argv[4:]]
        region(key, sys.argv[2], c, s - 1, e, [(a - 1, b) for a, b in marks]); return
    if len(sys.argv) > 1 and sys.argv[1] == "locus":
        key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
        os.makedirs(f"{R}/ucsc", exist_ok=True)
        locus(key, sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 10000); return
    if len(sys.argv) > 1 and sys.argv[1] == "overview":
        key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
        os.makedirs(f"{R}/ucsc", exist_ok=True)
        overview(key, f"{R}/ucsc/mafb_domain_overview.pdf", smooth=2)
        overview(key, f"{R}/ucsc/mafb_domain_overview_nosmooth.pdf", smooth=None); return
    win = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    first = int(sys.argv[2]) if len(sys.argv) > 2 else None
    key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
    S = pd.read_csv(f"{R}/mafb_domain/element_summaries.tsv", sep="\t")
    P = S.drop_duplicates(["chrom", "peak_start", "peak_end"]).sort_values("peak_start")[["cls", "chrom", "peak_start", "peak_end"]]
    if first: P = P.head(first)
    sc = os.environ.get("ACC_MAX", "30")
    out = f"{R}/ucsc/dar_peaks_win{win / 1000:g}k" + ("_pairscale" if sc == "pair" else f"_acc{sc}"); os.makedirs(out, exist_ok=True)
    bws = {n: pyBigWig.open(f"{W}/{n}-insertion-bin10-smooth20.bw") for n, _, _, _ in ACC}
    files = []
    for i, r in enumerate(P.itertuples(index=False), 1):
        mid = (r.peak_start + r.peak_end) // 2; a = max(0, mid - win // 2); b = a + win   # constant width
        lim = {}
        if os.environ.get("ACC_MAX", "30") != "pair":
            lim = {"RA": float(os.environ.get("ACC_MAX", "30")), "G4": float(os.environ.get("ACC_MAX", "30"))}
        else:
          for n, _, pair, _ in ACC:                       # pair scale: max over the view (bins at display resolution)
            m = np.nanmax(np.array(bws[n].stats(r.chrom, a, b, type="max", nBins=PIX), float))
            lim[pair] = max(lim.get(pair, 0), m)
          lim = {k: float(np.ceil(1.05 * v)) if v > 0 else 1.0 for k, v in lim.items()}
        cls = {"both": "shared"}.get(r.cls, r.cls)
        f = f"{out}/{i:02d}_{cls}_{r.chrom}_{r.peak_start + 1}-{r.peak_end}.pdf"
        render(key, r.chrom, a, b, r.peak_start, r.peak_end, lim, f); files.append(f)
        print(f"  {i:2d}/{len(P)} {cls:6s} {r.chrom}:{r.peak_start + 1:,}-{r.peak_end:,}  scale RA/Arco 0-{lim['RA']:g}, G4/G1 0-{lim['G4']:g}", flush=True)
        time.sleep(2)                                   # be gentle with the public server
    merged = f"{out}.pdf"
    subprocess.run(["pdfunite", *files, merged], check=True)
    print(f"  wrote {len(files)} PDFs to {out}/ and {merged}")


if __name__ == "__main__":
    main()
