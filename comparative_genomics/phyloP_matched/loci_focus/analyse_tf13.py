#!/usr/bin/env python
"""Paired motif-vs-own-enhancer constraint for the SCENIC+ MAFB regulators.

Same statistic as the AR analysis: each motif's phyloP minus the phyloP of the
remainder of the same DAR, obtained by subtraction so every motif is compared
against an identically-defined local background.
"""
import numpy as np, collections, json
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
rng=np.random.default_rng(0)
REG={"KLF12":"Glut-RA","KLF15":"Glut-RA","NR3C2":"Glut-RA","NFIA":"Glut-RA",
     "LIN28B":"Glut-RA","STAT1":"Glut-RA","FOSL2":"Glut-RA",
     "MAF":"GABA-4","SOX6":"GABA-4","TCF12":"GABA-4","ARX":"GABA-4",
     "SATB1":"GABA-4","MEIS2":"GABA-4"}
cls={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); cls[f[3]]=f[4]
mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
def load(tree):
    dar={}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p=l.split("\t"); dar[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    agg=collections.defaultdict(lambda:[0,0,0.0])
    for l in open(f"{F}/{tree}_tf13_hits.tab"):
        p=l.split("\t"); m,e,_=p[0].split("@")
        a=agg[(m,e)]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    out=collections.defaultdict(list)
    for (m,e),(msz,mcov,msum) in agg.items():
        if e not in dar or e not in cls: continue
        dsz,dcov,dsum=dar[e]
        rcov=dcov-mcov; rsum=dsum-msum
        if mcov<0.5*msz or rcov<50: continue
        loc="MAFB" if e in mafb else "other"
        out[(m,cls[e],loc)].append(msum/mcov-rsum/rcov)
    return out
def pv(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
def dod(a,b,n=20000):
    o=a.mean()-b.mean(); pool=np.concatenate([a,b]); k=len(a)
    p=np.array([(lambda q:q[:k].mean()-q[k:].mean())(rng.permutation(pool)) for _ in range(n)])
    return o,max((np.abs(p)>=abs(o)).mean(),1/n)
def get(D,m,c,loc=None):
    return np.array([x for (mm,cc,ll),v in D.items()
                     if mm==m and cc==c and (loc is None or ll==loc) for x in v])
AR={"subtree222":(0.044,0.021,0.023),"wholetree":(-0.033,-0.119,0.086)}
res={}
for tree,tl in (("subtree222","OSCINE SUBTREE"),("wholetree","WHOLE TREE")):
    D=load(tree)
    print(f"\n{'='*104}\n=== {tl}: motif minus its own enhancer, genome-wide ===")
    print(f"  {'TF':<8}{'SCENIC+':<10}{'Glut-RA':>10}{'p':>9}{'GABA-4':>10}{'p':>9}"
          f"{'DoD':>9}{'p':>9}{'n RA/G4':>14}")
    rows=[]
    for m in REG:
        a=get(D,m,"RA_over_Arco"); b=get(D,m,"GABA4_over_GABA1")
        if len(a)<40 or len(b)<40:
            print(f"  {m:<8}{REG[m]:<10}{'too few':>10}{'':>37}{str(len(a))+'/'+str(len(b)):>14}")
            continue
        o,pd_=dod(a,b)
        rows.append(dict(tf=m,reg=REG[m],ra=a.mean(),g4=b.mean(),dod=o,p=pd_,
                         pa=pv(a),pb=pv(b),na=len(a),nb=len(b)))
    for r in sorted(rows,key=lambda r:-r["dod"]):
        st=lambda p:"*" if p<.05 else " "
        print(f"  {r['tf']:<8}{r['reg']:<10}{r['ra']:>+10.3f}{format(r['pa'],'.4f')+st(r['pa']):>9}"
              f"{r['g4']:>+10.3f}{format(r['pb'],'.4f')+st(r['pb']):>9}"
              f"{r['dod']:>+9.3f}{format(r['p'],'.4f')+st(r['p']):>9}"
              f"{str(r['na'])+'/'+str(r['nb']):>14}")
    a,b,d=AR[tree]
    print(f"  {'AR':<8}{'(ref)':<10}{a:>+10.3f}{'':>9}{b:>+10.3f}{'':>9}{d:>+9.3f}")
    res[tree]=rows
    print(f"\n  --- at MAFB ---")
    print(f"  {'TF':<8}{'Glut-RA':>10}{'GABA-4':>10}{'DoD':>9}{'p':>9}{'n RA/G4':>12}")
    for m in REG:
        a=get(D,m,"RA_over_Arco","MAFB"); b=get(D,m,"GABA4_over_GABA1","MAFB")
        if len(a)<8 or len(b)<8:
            print(f"  {m:<8}{'too few':>10}{'':>28}{str(len(a))+'/'+str(len(b)):>12}"); continue
        o,p=dod(a,b)
        print(f"  {m:<8}{a.mean():>+10.3f}{b.mean():>+10.3f}{o:>+9.3f}"
              f"{format(p,'.4f')+('*' if p<.05 else ' '):>9}{str(len(a))+'/'+str(len(b)):>12}")
json.dump(res,open(f"{F}/results/compare_gw_mafb1p5mb/tf13_results.json","w"),indent=1,default=float)
