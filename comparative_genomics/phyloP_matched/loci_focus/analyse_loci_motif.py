#!/usr/bin/env python
"""Paired motif-vs-own-enhancer constraint at MAFB and EMX2, 192 archetypes.

Loci are the primary analysis; the genome-wide run provides the reference
distribution each locus value is placed against.
"""
import numpy as np, collections, json, os
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
rng=np.random.default_rng(0); MIN=12
meta={}
for l in open(f"{F}/loci_dars.bed"):
    f=l.rstrip("\n").split("\t"); meta[f[3]]=(f[4],f[5])
def load(tree,pre):
    dar={}
    for l in open(f"{F}/{tree}_{pre}_dars.tab"):
        p=l.split("\t"); dar[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    agg=collections.defaultdict(lambda:[0,0,0.0])
    for l in open(f"{F}/{tree}_{pre}_hits.tab"):
        p=l.split("\t"); m,e,_=p[0].split("@")
        a=agg[(m,e)]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    return dar,agg
def contrasts(dar,agg,keep=None):
    out=collections.defaultdict(list)
    for (m,e),(msz,mcov,msum) in agg.items():
        if e not in dar or e not in meta: continue
        cl,loc=meta[e]
        if keep and loc!=keep: continue
        dsz,dcov,dsum=dar[e]
        rcov=dcov-mcov; rsum=dsum-msum
        if mcov<0.5*msz or rcov<50: continue
        out[m].append((msum/mcov-rsum/rcov, cl))
    return out
def bp(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
def bh(ps):
    ps=np.array(ps); o=np.argsort(ps); q=ps[o]*len(ps)/np.arange(1,len(ps)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty(len(ps)); r[o]=np.clip(q,0,1); return r
def short(m):
    return m.split(".",2)[-1].replace("jaspar__","").replace("hocomoco__","")\
        .replace("swissregulon__mm__","").replace("swissregulon__hs__","")\
        .replace("cisbp__","").replace("transfac_public__","").replace("elemento__","")[:24]
gw={}
if os.path.exists(f"{F}/results/genome_wide/allmotif_subtree222.json"):
    for r in json.load(open(f"{F}/results/genome_wide/allmotif_subtree222.json")): gw[r["motif"]]=r
dar,agg=load("subtree222","lm")
for LOC in ("MAFB","EMX2"):
    C=contrasts(dar,agg,LOC)
    rows=[]
    for m,v in C.items():
        d=np.array([x[0] for x in v])
        if len(d)<MIN: continue
        ra=np.array([x[0] for x in v if x[1]=="RA_over_Arco"])
        g4=np.array([x[0] for x in v if x[1]=="GABA4_over_GABA1"])
        rows.append(dict(motif=m,n=len(d),d=d.mean(),p=bp(d),
                         d_ra=ra.mean() if len(ra)>=5 else np.nan,
                         d_g4=g4.mean() if len(g4)>=5 else np.nan,
                         n_ra=len(ra),n_g4=len(g4),
                         gw=gw.get(m,{}).get("d_both",np.nan)))
    q=bh([r["p"] for r in rows])
    for r,v in zip(rows,q): r["q"]=v
    rows.sort(key=lambda r:-r["d"])
    print(f"\n{'='*104}\n=== {LOC}: motif vs the rest of its own enhancer, oscine subtree "
          f"({len(rows)} motifs with >={MIN} elements) ===")
    print(f"  {'rank':>5}  {'motif archetype':<26}{'n':>5}{'diff':>9}{'q':>8}"
          f"{'Glut-RA':>10}{'GABA-4':>9}{'genome-wide':>13}")
    AR=[i for i,r in enumerate(rows) if "MA0007.3" in r["motif"]]
    for i,r in enumerate(rows):
        if i<10 or i>=len(rows)-4 or i in AR:
            t="  <- AR" if i in AR else ""
            f_=lambda v:f"{v:+.3f}" if np.isfinite(v) else "  --  "
            print(f"  {i+1:>5}  {short(r['motif']):<26}{r['n']:>5}{r['d']:>+9.3f}{r['q']:>8.3f}"
                  f"{f_(r['d_ra']):>10}{f_(r['d_g4']):>9}{f_(r['gw']):>13}{t}")
        elif i==10: print(f"  {'...':>5}")
    s=[r for r in rows if r["q"]<0.05]
    print(f"  FDR 5%: {len(s)}/{len(rows)} motifs   (positive {sum(1 for r in s if r['d']>0)})")
    if AR:
        r=rows[AR[0]]
        print(f"  AR rank {AR[0]+1}/{len(rows)}  d={r['d']:+.3f} q={r['q']:.3f}"
              f"  | genome-wide d={r['gw']:+.3f}" if np.isfinite(r["gw"]) else "")
    json.dump(rows,open(result(f"locimotif_{LOC}.json", LOC),"w"),indent=1,default=float)
