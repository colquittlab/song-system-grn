#!/usr/bin/env python
"""Single source for the AR report: regenerates every table from the saved data."""
import json, pickle, collections, sys, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"maf_switch")
import ar_core as C
from ar_highic import IC, HI, pairs_ic
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
GS,GE=13365864,13369534; DS,DE=12551036,13547533
rng=np.random.default_rng(0)
def p1(a,n=20000):
    s=rng.choice([-1,1],(n,len(a)))
    return float(a.mean()),float(((np.abs((s*a).mean(1))>=abs(a.mean())).sum()+1)/(n+1))
def perm(a,b,n=20000):
    o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0
    for _ in range(n):
        p=rng.permutation(c); h+=abs(p[:k].mean()-p[k:].mean())>=abs(o)
    return float(o),float((h+1)/(n+1))
cls={}
for l in open("allelem_cls.bed"):
    g=l.rstrip("\n").split("\t"); cls[g[3]]=(g[0],int(g[1]),int(g[2]),g[4])
m={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
d=lambda e:(lambda ch,s,t,c: 0 if (t>GS and s<GE) else (s-GE if s>=GE else GS-t))(*cls[e])
SETS={"domain":{e for e in m if cls[e][0]=="chr20" and cls[e][2]>DS and cls[e][1]<DE},
      "500kb":{e for e in m if d(e)<=500_000},
      "1Mb":{e for e in m if d(e)<=1_000_000},
      "1.5Mb":set(m)}
OUT={"sets":{k:{"n":len(v),
                "ra":sum(1 for e in v if cls[e][3]=="RA_over_Arco"),
                "g4":sum(1 for e in v if cls[e][3]!="RA_over_Arco")} for k,v in SETS.items()}}
def cell(P,tag):
    a=C.group(P,"ra",tag); b=C.group(P,"g4",tag)
    if len(a)<3 or len(b)<3: return None
    av=a[:,0]-a[:,1]; bv=b[:,0]-b[:,1]
    ra,pra=p1(av); g4,pg4=p1(bv); dd,pd=perm(av,bv)
    return dict(ra=ra,p_ra=pra,g4=g4,p_g4=pg4,dod=dd,p_dod=pd,n_ra=len(av),n_g4=len(bv))
# --- part 1a: genome-wide, both trees, all/high/low IC
OUT["gw"]={}
for tree in ("subtree222","wholetree"):
    for cols in ("all","hi","lo"):
        for mk in ("full","half","jaspar"):
            P=pairs_ic(tree,mk,keep=None,cols=cols)
            r=cell(P,None)
            if r: OUT["gw"][f"{tree}|{cols}|{mk}"]=r
# --- part 1b: distance series, oscine, all columns
OUT["dist"]={}
for lab,S in SETS.items():
    for mk in ("full","half","jaspar"):
        r=cell(C.pairs("subtree222",mk,"site",keep=S),"sub")
        if r: OUT["dist"][f"{lab}|{mk}"]=r
# --- part 1c: domain, column-restricted
OUT["dom_ic"]={}
for cols in ("all","hi","lo"):
    for mk in ("full","half","jaspar"):
        P=pairs_ic("subtree222",mk,keep=SETS["domain"],cols=cols)
        r=cell(P,"sub")
        if r: OUT["dom_ic"][f"{cols}|{mk}"]=r
# --- part 1d: motif density per Mb of DAR
bp=collections.Counter()
for l in open("subtree222_am_dars.tab"):
    p=l.split("\t")
    if p[0] in cls: bp[cls[p[0]][3]]+=int(p[2])
site=collections.Counter()
for l in open("ar3_sites.bed"):
    mo,s,loc,e,_=l.rstrip("\n").split("\t")[3].split("|")
    if e in cls: site[(mo,cls[e][3])]+=1
OUT["density"]={mo:{"ra":site[(mo,"RA_over_Arco")]/(bp["RA_over_Arco"]/1e6),
                    "g4":site[(mo,"GABA4_over_GABA1")]/(bp["GABA4_over_GABA1"]/1e6)}
                for mo in ("full","half","jaspar")}
OUT["ic"]={k:[float(x) for x in v] for k,v in IC.items()}
json.dump(OUT,open("results/compare_mafb_spans/ar_report_part1.json","w"),indent=1)
print("  part 1 written:",len(OUT["gw"]),"gw cells,",len(OUT["dist"]),"distance cells,",
      len(OUT["dom_ic"]),"domain cells")
