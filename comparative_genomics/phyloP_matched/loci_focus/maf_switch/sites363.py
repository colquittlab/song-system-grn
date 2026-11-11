#!/usr/bin/env python
"""Per-site, per-column base composition from the avian363 (lonStrDom2-anchored) MAF.

Same alignment the phyloP tracks were computed on, so the oscine/outgroup split
here is the same clade the subtree222 track uses: birdAnc222 = Passeri,
146 oscines vs 217 outgroup.

Writes  data/cols363.json   {refpos: {genome: base}}  (ancestors included)
        data/sites363.json  per site, per motif column, composition summary
"""
import collections, json, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
D=f"{F}/maf_switch"
COMP=str.maketrans("ACGTacgtNn-","TGCAtgcaNn-")
OSC={l.strip() for l in open(f"{D}/data/osc363.txt") if l.strip()}
OUT={l.strip() for l in open(f"{D}/data/out363.txt") if l.strip()}
REF="Lonchura_striata"

def load(path):
    col=collections.defaultdict(dict); block=[]; tree=None
    def flush(b):
        ref=None
        for sp,st,txt in b:
            if sp.split(".")[0]==REF: ref=(st,txt); break
        if ref is None: return
        st,rtxt=ref; p=st; idx=[]
        for i,ch in enumerate(rtxt):
            if ch!="-": idx.append((i,p)); p+=1
        for sp,s2,txt in b:
            g=sp.split(".")[0]
            for i,rp in idx:
                if i<len(txt):
                    c=txt[i].upper()
                    col[rp].setdefault(g,c if c in "ACGT-" else "N")
    for l in open(path):
        if l.startswith("# hal"): tree=l[6:].strip()
        elif l.startswith("a"): flush(block); block=[]
        elif l.startswith("s\t") or l.startswith("s "):
            f=l.rstrip("\n").split("\t") if "\t" in l else l.split()
            block.append((f[1],int(f[2]),f[6]))
    flush(block)
    return col,tree

def summarise(d):
    def maj(sp):
        c=collections.Counter(b for s,b in d.items() if s in sp and b in "ACGT")
        n=sum(c.values())
        if not n: return ("-",0.0,0,{})
        b,k=c.most_common(1)[0]; return (b,k/n,n,dict(c))
    return maj(OSC),maj(OUT)

if __name__=="__main__":
    maf=sys.argv[1] if len(sys.argv)>1 else f"{D}/data/ar_mafb_363.maf"
    col,tree=load(maf)
    print(f"  reference columns: {len(col):,}   genomes seen: "
          f"{len({g for d in col.values() for g in d}):,}")
    if tree: open(f"{D}/data/tree363.nw","w").write(tree+"\n")
    json.dump({str(k):v for k,v in col.items()},open(f"{D}/data/cols363.json","w"))
    beds=[(f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","full",12),
          (f"{F}/tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed","jaspar",16)]
    rec=[]
    for bed,tag,W in beds:
        for l in open(bed):
            f=l.rstrip("\n").split("\t")
            start,end,name,strand=int(f[1]),int(f[2]),f[3],f[5]
            cols=[]
            for k in range(W):
                rp=start+k if strand=="+" else end-1-k
                d=col.get(rp,{})
                if strand=="-": d={s:b.translate(COMP) for s,b in d.items()}
                (om,of,on,oc),(gm,gf,gn,gc)=summarise(d)
                cols.append(dict(k=k+1,refpos=rp,osc=om,oscF=round(of,3),oscN=on,oscC=oc,
                                 out=gm,outF=round(gf,3),outN=gn,outC=gc))
            rec.append(dict(motif=tag,site=name,strand=strand,start=start,end=end,cols=cols))
    json.dump(rec,open(f"{D}/data/sites363.json","w"),indent=1)
    print(f"  sites summarised: {len(rec)}")
