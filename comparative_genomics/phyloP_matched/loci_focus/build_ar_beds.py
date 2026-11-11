#!/usr/bin/env python
"""Hit and merged-site BEDs for the three AR motifs, in one naming scheme.

name = motif|class|locus|element|index
Hits are converted from element-relative FIMO coordinates to genome coordinates.
Merged sites collapse overlapping hits within an element, which for a
palindromic motif such as the ARE means collapsing the two strand reports of
one site into one.
"""
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
MOT={"full":"motif2","half":"motif4","jaspar":"jaspar"}
el={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); el[f[3]]=(f[0],int(f[1]),int(f[2]),f[4])
mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
hb=open(f"{F}/ar3_hits.bed","w"); mb=open(f"{F}/ar3_sites.bed","w")
nh=collections.Counter(); ns=collections.Counter()
for tag,pref in MOT.items():
    per=collections.defaultdict(list)
    for s in ("ra","g4"):
        for i,l in enumerate(open(f"{F}/fimo_{pref}_{s}.txt")):
            if not i: continue
            p=l.rstrip("\n").split("\t"); e=p[2]
            if e not in el: continue
            ch,es,ee,cl=el[e]
            if (cl=="RA_over_Arco")!=(s=="ra"): continue
            per[(s,e)].append((es+int(p[3])-1, es+int(p[4])))
    for (s,e),v in per.items():
        ch,es,ee,cl=el[e]; loc="MAFB" if e in mafb else "other"
        for j,(a,b) in enumerate(sorted(v)):
            hb.write(f"{ch}\t{a}\t{b}\t{tag}|{s}|{loc}|{e}|{j}\n"); nh[tag]+=1
        mg=[]
        for a,b in sorted(v):
            if mg and a<mg[-1][1]: mg[-1][1]=max(mg[-1][1],b)
            else: mg.append([a,b])
        for j,(a,b) in enumerate(mg):
            mb.write(f"{ch}\t{a}\t{b}\t{tag}|{s}|{loc}|{e}|{j}\n"); ns[tag]+=1
hb.close(); mb.close()
for t in MOT: print(f"  {t:7s} hits {nh[t]:>7,}   merged sites {ns[t]:>7,}")
