#!/usr/bin/env python
"""Per-position BEDs for the three AR motifs, one representative hit per site.

A palindromic motif is reported by FIMO on both strands at the same locus. Both
reports carry strand-corrected position indices, so keeping both would fold
column k onto column W-1-k and symmetrise the profile. Overlapping hits are
therefore merged and only the best-scoring hit of each site is emitted.

name = motif|class|locus|element|site|k
"""
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
MOT={"full":("motif2",12),"half":("motif4",8),"jaspar":("jaspar",16)}
el={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); el[f[3]]=(f[0],int(f[1]),int(f[2]),f[4])
mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
out=open(f"{F}/ar3_pos.bed","w"); kept=collections.Counter(); drop=collections.Counter()
for tag,(pref,W) in MOT.items():
    per=collections.defaultdict(list)
    for s in ("ra","g4"):
        for i,l in enumerate(open(f"{F}/fimo_{pref}_{s}.txt")):
            if not i: continue
            p=l.rstrip("\n").split("\t"); e=p[2]
            if e not in el: continue
            ch,es,ee,cl=el[e]
            if (cl=="RA_over_Arco")!=(s=="ra"): continue
            st,sp,strand,sc=int(p[3]),int(p[4]),p[5],float(p[6])
            if sp-st+1!=W: continue                   # guard against ragged hits
            per[(s,e)].append((es+st-1,es+sp,strand,sc))
    for (s,e),v in per.items():
        ch,_,_,_=el[e]; loc="MAFB" if e in mafb else "other"
        v.sort()
        groups=[]
        for h in v:
            if groups and h[0]<groups[-1][-1][1]: groups[-1].append(h)
            else: groups.append([h])
        for j,g in enumerate(groups):
            a,b,strand,sc=max(g,key=lambda x:x[3])    # best-scoring hit of the site
            drop[tag]+=len(g)-1; kept[tag]+=1
            for k in range(W):
                pos = a+k if strand=="+" else b-1-k
                out.write(f"{ch}\t{pos}\t{pos+1}\t{tag}|{s}|{loc}|{e}|{j}|{k}\n")
out.close()
for t in MOT: print(f"  {t:7s} sites kept {kept[t]:>7,}   duplicate hits dropped {drop[t]:>7,}")
