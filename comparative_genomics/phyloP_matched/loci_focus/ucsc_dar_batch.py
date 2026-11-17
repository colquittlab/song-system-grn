#!/usr/bin/env python
"""UCSC Genome Browser PDFs of each MAFB-domain DAR peak, from the lab hub's own bigWigs.

For every SnapATAC2 peak in results/mafb_domain/element_summaries.tsv (RA, G4 and shared elements;
one view per peak, lobes not split), the browser draws the peak +/- FLANK with only these tracks, in
this order, as custom tracks pointing at the hub's bigWig files (so order, colours and scales are set
here rather than by the hub's trackDb):
  Glut-RA, Glut-Arco-1, GABA-4-1, GABA-1-1 accessibility (insertion, 10-bp bins, 20-bp smoothing)
  oscine-lineage phyloP (subtree birdAnc222) and whole-tree phyloP (363 birds), per base
Glut-RA/Glut-Arco-1 share a fixed scale, as do GABA-4-1/GABA-1-1: 0 to 1.05 x the pair's maximum over
the view. phyloP limits follow the hub (oscine -3:3, whole tree -4:6). Vertical guidelines and the
descriptions above tracks are off; the peak is highlighted in grey. The UCSC API key is read from
~/.ucsc_api_key (never stored here).

    python ucsc_dar_batch.py [flank_bp=10000] [first N peaks]
Writes results/ucsc/dar_peaks_flank<F>/<nn>_<class>_chr20_<start>-<end>.pdf and dar_peaks_flank<F>.pdf
(all pages, in genome order).
"""
import os, re, subprocess, sys, time, urllib.parse, urllib.request, http.cookiejar
import numpy as np, pandas as pd, pyBigWig

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results"
BASE = "https://colquitt-track-hub.pbsci.ucsc.edu/lonStrDom2"
HUB, GENOME, UCSC = f"{BASE}/hub.txt", "GCF_005870125.1", "https://genome.ucsc.edu"
W = f"{BASE}/{GENOME}/scATAC/motor-pathway/wigs/insertion_bin10-smooth20"
ACC = [("Glut-RA", "227,26,28", "RA"), ("Glut-Arco-1", "251,154,153", "RA"),
       ("GABA-4-1", "0,0,0", "G4"), ("GABA-1-1", "0,0,0", "G4")]
PHY = [("phyloP oscine", f"{BASE}/{GENOME}/comparative/subtree222_perbase_lsd2.bw", "63,107,69", "154,123,36", "-3:3"),
       ("phyloP all 363", f"{BASE}/{GENOME}/comparative/wholetree_perbase_lsd2.bw", "33,102,172", "178,24,43", "-4:6")]
PIX = 1200


def custom_text(lim):
    t = []
    for name, col, pair in ACC:
        t.append(f'track type=bigWig name="{name} ATAC" bigDataUrl={W}/{name}-insertion-bin10-smooth20.bw color={col} '
                 f'autoScale=off viewLimits=0:{lim[pair]:g} alwaysZero=on maxHeightPixels=128:30:11 visibility=full')
    for name, url, col, alt, vl in PHY:
        t.append(f'track type=bigWig name="{name}" bigDataUrl={url} color={col} altColor={alt} autoScale=off '
                 f'viewLimits={vl} alwaysZero=on windowingFunction=mean maxHeightPixels=128:30:11 visibility=full')
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
    flank = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    first = int(sys.argv[2]) if len(sys.argv) > 2 else None
    key = open(os.path.expanduser("~/.ucsc_api_key")).read().strip()
    S = pd.read_csv(f"{R}/mafb_domain/element_summaries.tsv", sep="\t")
    P = S.drop_duplicates(["chrom", "peak_start", "peak_end"]).sort_values("peak_start")[["cls", "chrom", "peak_start", "peak_end"]]
    if first: P = P.head(first)
    out = f"{R}/ucsc/dar_peaks_flank{flank // 1000}k"; os.makedirs(out, exist_ok=True)
    bws = {n: pyBigWig.open(f"{W}/{n}-insertion-bin10-smooth20.bw") for n, _, _ in ACC}
    files = []
    for i, r in enumerate(P.itertuples(index=False), 1):
        a, b = max(0, r.peak_start - flank), r.peak_end + flank
        lim = {}
        for n, _, pair in ACC:                       # pair scale: max over the view (bins at display resolution)
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
