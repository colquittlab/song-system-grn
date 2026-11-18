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


def custom_text(lim):
    t = []
    for name, col, pair, lab in ACC:
        t.append(f'track type=bigWig name="{lab}" bigDataUrl={W}/{name}-insertion-bin10-smooth20.bw color={col} '
                 f'autoScale=off viewLimits=0:{lim[pair]:g} alwaysZero=on maxHeightPixels=128:{HEIGHT}:11 visibility=full')
    for name, url, col, alt, vl in PHY:
        t.append(f'track type=bigWig name="{name}" bigDataUrl={url} color={col} altColor={alt} autoScale=off '
                 f'viewLimits={vl} alwaysZero=on windowingFunction=mean maxHeightPixels=128:{HEIGHT}:11 visibility=full')
    return "\n".join(t)


def render(key, chrom, a, b, hs, he, lim, out):
    """two requests in one session: load the custom tracks (their ids are assigned by the server), then
    hide every track and show only those ids -- hideTracks would otherwise hide the custom tracks too"""
    cj = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    base = {"hubUrl": HUB, "genome": GENOME, "position": f"{chrom}:{a + 1}-{b}", "pix": PIX, "apiKey": key}
    post = lambda q: cj.open(f"{UCSC}/cgi-bin/hgTracks", urllib.parse.urlencode(q).encode(), timeout=300).read().decode("utf-8", "replace")
    for i in range(3):
        try:
            page = post({**base, "hgt.customText": custom_text(lim)})
            ids = list(dict.fromkeys(re.findall(r"ct_[A-Za-z0-9]+_[0-9]+", page)))
            if len(ids) < len(ACC) + len(PHY): raise RuntimeError(f"custom tracks not loaded ({ids})")
            page = post({**base, "hideTracks": 1, **{t: "full" for t in ids}, "guidelines": "off", "centerLabels": "off",
                         "highlight": f"{chrom}:{hs + 1}-{he}#e6e6e6", "hgt.psOutput": "on"})
            p = re.search(r"\.\./trash/hgt/hgt_genome_[^'\"]*\.pdf", page)
            if not p: raise RuntimeError("no PDF link in response")
            open(out, "wb").write(cj.open(f"{UCSC}/{p.group(0)[3:]}", timeout=120).read()); return
        except Exception:
            if i == 2: raise
            time.sleep(10 * (i + 1))


def main():
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
