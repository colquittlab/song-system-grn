#!/usr/bin/env python
"""Per-site, per-column base composition across species from the lonStrDom2 MAF.

For each motif site, walks the reference row of every overlapping MAF block to
assign reference coordinates to alignment columns, then reads each species'
base in that column. Motif position is strand-corrected and complemented for
minus-strand hits, so column 1 is always the motif's first column.

Ancestral (Anc*) rows are skipped; only extant genomes are counted.
"""
import collections, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
COMP=str.maketrans("ACGTacgtNn-","TGCAtgcaNn-")
OSC={l.strip() for l in open(f"{F}/maf_oscines.txt") if l.strip()}
OUT={l.strip() for l in open(f"{F}/maf_outgroup.txt") if l.strip()}

def load_maf(path):
    """-> {refpos: {species: base}}"""
    col=collections.defaultdict(dict)
    block=[]
    def flush(b):
        if not b: return
        ref=None
        for sp,st,sz,strand,txt in b:
            if sp.startswith("Lonchura_striata"): ref=(st,txt); break
        if ref is None: return
        st,rtxt=ref
        idx=[]; p=st
        for i,ch in enumerate(rtxt):
            if ch!="-": idx.append((i,p)); p+=1
        for sp,s2,z2,strand,txt in b:
            g=sp.split(".")[0]
            if g.startswith("Anc"): continue
            for i,rp in idx:
                if i<len(txt):
                    c=txt[i].upper()
                    if c in "ACGT": col[rp].setdefault(g,c)
                    elif c=="-": col[rp].setdefault(g,"-")
    for l in open(path):
        if l.startswith("a "): flush(block); block=[]
        elif l.startswith("s "):
            f=l.split()
            block.append((f[1],int(f[2]),int(f[3]),f[4],f[6]))
    flush(block)
    return col

def site_table(col,start,end,strand,W):
    """-> list over motif positions of {species: base}"""
    out=[]
    for k in range(W):
        rp = start+k if strand=="+" else end-1-k
        d=col.get(rp,{})
        if strand=="-": d={s:b.translate(COMP) for s,b in d.items()}
        out.append(d)
    return out

def summarise(d):
    """-> (osc_major, osc_freq, osc_n, out_major, out_freq, out_n)"""
    def maj(sp):
        c=collections.Counter(b for s,b in d.items() if s in sp and b in "ACGT")
        n=sum(c.values())
        if not n: return ("-",0.0,0)
        b,k=c.most_common(1)[0]; return (b,k/n,n)
    return maj(OSC)+maj(OUT)

if __name__=="__main__":
    col=load_maf(sys.argv[1] if len(sys.argv)>1 else "/tmp/ar_mafb.maf")
    print(f"  reference columns loaded: {len(col):,}")
    rows=[]
    for l in open("/tmp/ar_mafb.bed"):
        f=l.rstrip("\n").split("\t")
        start,end=int(f[1]),int(f[2]); tag,name,strand=f[3].split("|")
        W=end-start
        tab=site_table(col,start,end,strand,W)
        rows.append((tag,name,strand,start,end,tab))
    import json
    rec=[]
    for tag,name,strand,start,end,tab in rows:
        cols=[]
        for k,d in enumerate(tab):
            om,of,on,gm,gf,gn=summarise(d)
            co={b:sum(1 for s_,x in d.items() if s_ in OSC and x==b) for b in "ACGT"}
            cg={b:sum(1 for s_,x in d.items() if s_ in OUT and x==b) for b in "ACGT"}
            cols.append(dict(k=k+1,osc=om,oscF=round(of,3),oscN=on,
                             out=gm,outF=round(gf,3),outN=gn,
                             oscC={b:n for b,n in co.items() if n},
                             outC={b:n for b,n in cg.items() if n}))
        rec.append(dict(motif=tag,site=name,strand=strand,start=start,end=end,cols=cols))
    json.dump(rec,open(f"{F}/results/mafb_1p5mb/maf_ar_mafb.json","w"),indent=1)
    print(f"  sites summarised: {len(rec)}")
