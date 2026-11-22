#!/usr/bin/env python
"""Is the Glut-RA / GABA-4 switch asymmetry more than chance?

Null: a site's DAR class is exchangeable with respect to how many of its
columns switch. Sites differ in how many columns are testable (alignment
coverage), so the statistic is a POOLED RATE -- switches per testable column --
not a per-site mean, and the permutation reshuffles class labels across sites
while carrying each site's testable count with it.
"""
import json, sys, numpy as np, collections
NPERM=20000
def load(path,minO,minG,thr=0.80):
    R=json.load(open(path)); out=collections.defaultdict(list)
    for r in R:
        cls="Glut-RA" if r["site"].startswith("RA_") else "GABA-4"
        test=sw=0
        for c in r["cols"]:
            if c["oscN"]<minO or c["outN"]<minG: continue
            test+=1
            if c["osc"]!=c["out"] and c["oscF"]>=thr and c["outF"]>=thr: sw+=1
        out[r["motif"]].append((cls,sw,test))
    return out
def test(rows,rng):
    cls=np.array([r[0]=="Glut-RA" for r in rows])
    sw=np.array([r[1] for r in rows],float); te=np.array([r[2] for r in rows],float)
    def stat(m):
        a=te[m].sum(); b=te[~m].sum()
        if a==0 or b==0: return np.nan
        return sw[m].sum()/a - sw[~m].sum()/b
    obs=stat(cls); n=len(rows); k=cls.sum(); hit=0
    for _ in range(NPERM):
        p=np.zeros(n,bool); p[rng.choice(n,k,replace=False)]=True
        if abs(stat(p))>=abs(obs)-1e-12: hit+=1
    return obs,(hit+1)/(NPERM+1),sw[cls].sum(),te[cls].sum(),sw[~cls].sum(),te[~cls].sum()
for path,lab,minO,minG in ((sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])),):
    D=load(path,minO,minG)
    print(f"\n=== {lab}   (coverage floor: {minO} oscine / {minG} outgroup)")
    print(f"  {'motif':<10}{'RA switches/cols':>20}{'G4 switches/cols':>20}"
          f"{'rate diff':>12}{'p':>9}")
    for mk in ("full","jaspar"):
        rows=D.get(mk,[])
        if not rows: continue
        rng=np.random.default_rng(0)
        obs,p,sa,ta,sb,tb=test(rows,rng)
        print(f"  {mk:<10}{f'{sa:.0f}/{ta:.0f}':>20}{f'{sb:.0f}/{tb:.0f}':>20}"
              f"{obs:>+12.4f}{p:>9.4f}{'*' if p<.05 else ''}")
    # motifs pooled, for context only -- the two site sets overlap in sequence
    rows=[r for mk in ("full","jaspar") for r in D.get(mk,[])]
    rng=np.random.default_rng(0)
    obs,p,sa,ta,sb,tb=test(rows,rng)
    print(f"  {'pooled':<10}{f'{sa:.0f}/{ta:.0f}':>20}{f'{sb:.0f}/{tb:.0f}':>20}"
          f"{obs:>+12.4f}{p:>9.4f}{'*' if p<.05 else ''}   (motifs overlap; context only)")
