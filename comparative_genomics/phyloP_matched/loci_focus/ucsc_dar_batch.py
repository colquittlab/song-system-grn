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
HEIGHT = 45                                          # track height (px): 1.5 x the hub default of 30


GENES = f"{BASE}/{GENOME}/genes/GCF_005870125.1_lonStrDom2_genomic_stringtie_FPKMthresh0.5_FPKMmaxfrac0.1_minhits5_slim.bb"


def gene_track(chrom, a, b):
    """inline BED12 of the hub's gene models over [a, b): the longest transcript per gene, named by symbol (name2)"""
    bb = pyBigWig.open(GENES); rows = {}
    for s, e, rest in bb.entries("NC_042585.1" if chrom == "chr20" else chrom, a, b):
        f = rest.split("\t"); sym = f[9] or f[0]
        if sym not in rows or e - s > rows[sym][1] - rows[sym][0]: rows[sym] = (s, e, f)
    out = ['track name="genes" description="NCBI gene models (3\' UTR extended), longest transcript per gene" '
           'visibility=pack color=60,60,60']
    for sym, (s, e, f) in sorted(rows.items(), key=lambda kv: kv[1][0]):
        out.append("\t".join(map(str, [chrom, s, e, sym, 0, f[2], f[3], f[4], 0, f[6], f[7], f[8]])))
    return "\n".join(out)


NATIVE = {"chr20": "NC_042585.1"}
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


def custom_text(lim, phy=True, height=None, auto=False, extra="", pre=""):
    height = height or HEIGHT
    t = []
    for name, col, pair, lab in ACC:
        t.append(f'track type=bigWig name="{lab}" bigDataUrl={W}/{name}-insertion-bin10-smooth20.bw color={col} '
                 + ("autoScale=on " if auto else f"autoScale=off viewLimits=0:{lim[pair]:g} ")
                 + f'alwaysZero=on maxHeightPixels=128:{height}:11 visibility=full')
    for name, url, col, alt, vl in (PHY if phy else []):
        t.append(f'track type=bigWig name="{name}" bigDataUrl={url} color={col} altColor={alt} autoScale=off '
                 f'viewLimits={vl} alwaysZero=on windowingFunction=mean maxHeightPixels=128:{HEIGHT}:11 visibility=full')
    return (pre + "\n" if pre else "") + "\n".join(t) + ("\n" + extra if extra else "")


def render(key, chrom, a, b, hs, he, lim, out, phy=True, highlight=None, height=None, auto=False, show=None, extra="", pre=""):
    """show: hub tracks to keep visible as {track: visibility}; highlight "" = none"""
    """two requests in one session: load the custom tracks (their ids are assigned by the server), then
    hide every track and show only those ids -- hideTracks would otherwise hide the custom tracks too"""
    cj = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    base = {"hubUrl": HUB, "genome": GENOME, "position": f"{chrom}:{a + 1}-{b}", "pix": PIX, "apiKey": key}
    post = lambda q: cj.open(f"{UCSC}/cgi-bin/hgTracks", urllib.parse.urlencode(q).encode(), timeout=300).read().decode("utf-8", "replace")
    for i in range(3):
        try:
            page = post({**base, "hgt.customText": custom_text(lim, phy, height, auto, extra, pre)})
            ids = list(dict.fromkeys(re.findall(r"ct_[A-Za-z0-9]+_[0-9]+", page)))
            if len(ids) < len(ACC) + (len(PHY) if phy else 0) + bool(extra) + bool(pre): raise RuntimeError(f"custom tracks not loaded ({ids})")
            q = {**base, "hideTracks": 1, **{t: ("pack" if t.startswith("ct_genes_") else "full") for t in ids}, **(show or {}),
                 "guidelines": "off", "centerLabels": "off",
                 "hgt.psOutput": "on"}
            if highlight != "": q["highlight"] = highlight or f"{chrom}:{hs + 1}-{he}#e6e6e6"
            page = post(q)
            p = re.search(r"\.\./trash/hgt/hgt_genome_[^'\"]*\.pdf", page)
            if not p: raise RuntimeError("no PDF link in response")
            open(out, "wb").write(cj.open(f"{UCSC}/{p.group(0)[3:]}", timeout=120).read()); return
        except Exception:
            if i == 2: raise
            time.sleep(10 * (i + 1))


def overview(key, out):
    """the whole DHX35-TOP1 domain: RA + C1H-1 peak-to-gene links to MAFB (1-Mb run, r >= 0.45) above the
    accessibility tracks (autoscaled, 30 px), then the hub's gene models"""
    pre, n = p2g_track("chr20", 12551035, 13547533)
    render(key, "chr20", 12551035, 13547533, None, None, {}, out, phy=False, highlight="", height=30, auto=True,
           extra=gene_track("chr20", 12551035, 13547533), pre=pre)
    print(f"  RA + C1H-1 links to MAFB drawn: {n}")
    print(f"  wrote {out}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "overview":
        key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
        os.makedirs(f"{R}/ucsc", exist_ok=True)
        overview(key, f"{R}/ucsc/mafb_domain_overview.pdf"); return
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
