#!/usr/bin/env python
"""Per-base alignment depth by clade over a reference region, from the 363-way MAF chunks.

For every Lonchura_striata reference base in the region: the fraction of oscine (Anc222 clade, 145
species, reference excluded) and of non-oscine (217) species with a non-gap base aligned to it.
A species with several rows in a block (paralogous copies) counts once. Reference bases in no block
get 0. Clade lists are those of the clade phastCons run (phastcons_clade/{osc,nonosc}.txt).

    python domain_depth.py chr20 12541036 13557533 chunk.maf [chunk.maf ...]
Writes results/mafb_domain/depth_{osc,nonosc}.bw (whole chromosome length; values only in the region).
"""
import sys
import numpy as np, pyBigWig

P = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched"
OUT = f"{P}/loci_focus/results/mafb_domain"
REF = "Lonchura_striata"


def main():
    chrom, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); mafs = sys.argv[4:]
    osc = {l.strip() for l in open(f"{P}/phastcons_clade/osc.txt") if l.strip()} - {REF}
    non = {l.strip() for l in open(f"{P}/phastcons_clade/nonosc.txt") if l.strip()}
    cnt = {"osc": np.zeros(b - a, np.int16), "nonosc": np.zeros(b - a, np.int16)}
    size = None
    for f in mafs:
        blk = []
        def flush():
            nonlocal size
            ref = [x for x in blk if x[0] == REF]
            if not ref: return
            _, s, strand, src, seq = ref[0]; size = src
            if strand != "+" or s >= b: return
            cols = np.frombuffer(seq.encode(), "S1") != b"-"; pos = s + np.cumsum(cols) - 1   # ref coordinate per column
            keep = cols & (pos >= a) & (pos < b)
            if not keep.any(): return
            have = {}                                  # species -> columns with a base in any of its rows
            for name, _, _, _, sq in blk:
                if name in osc or name in non:
                    m = np.frombuffer(sq.encode(), "S1") != b"-"
                    have[name] = have[name] | m if name in have else m
            for clade, sp in (("osc", osc), ("nonosc", non)):
                hits = np.zeros(len(seq), np.int16)
                for name, m in have.items():
                    if name in sp: hits += m
                cnt[clade][pos[keep] - a] += hits[keep]
        for line in open(f):
            if line.startswith("a"):
                if blk: flush()
                blk = []
            elif line.startswith("s "):
                x = line.split(); sp, ctg = x[1].split(".", 1)
                if sp == REF and ctg != chrom: blk = [("__skip__", 0, "+", 0, "")]; continue
                if blk and blk[0][0] == "__skip__": continue
                blk.append((sp, int(x[2]), x[4], int(x[5]), x[6]))
        if blk and blk[0][0] != "__skip__": flush()
        print(f"  {f}: done", flush=True)
    for clade, n in (("osc", len(osc)), ("nonosc", len(non))):
        v = cnt[clade] / n
        bw = pyBigWig.open(f"{OUT}/depth_{clade}.bw", "w"); bw.addHeader([(chrom, size)])
        bw.addEntries(chrom, a, values=v.astype(float).tolist(), span=1, step=1); bw.close()
        print(f"  depth_{clade}.bw: mean {v.mean():.3f}, bases < 0.5: {(v < 0.5).mean():.1%}")


if __name__ == "__main__":
    main()
