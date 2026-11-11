#!/usr/bin/env python
"""AR constraint at MAFB as a function of distance from the gene.

The published "MAFB locus" set is the gene +/- 1.5 Mb (a +/-500 kb window then
expanded by 1 Mb for power). This asks whether the result holds when the set is
pulled back toward the gene.
"""
import sys, pickle, numpy as np
sys.path.insert(0,".")
import ar_core as C
rng=np.random.default_rng(0)
def p1(a,n=20000):
    s=rng.choice([-1,1],(n,len(a)))
    return a.mean(),((np.abs((s*a).mean(1))>=abs(a.mean())).sum()+1)/(n+1)
def perm(a,b,n=20000):
    o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0
    for _ in range(n):
        p=rng.permutation(c); h+=abs(p[:k].mean()-p[k:].mean())>=abs(o)
    return o,(h+1)/(n+1)
SETS=pickle.load(open("/tmp/mafb_sets.pkl","rb"))
st=lambda p:"*" if p<.05 else " "
for tree,tl in (("subtree222","oscine"),("wholetree","whole tree")):
    print(f"\n=== {tl}: motif minus rest of its own enhancer, merged-site unit")
    print(f"  {'set':<16}{'motif':<9}{'Glut-RA':>19}{'GABA-4':>19}{'DoD':>19}   n")
    for lab,S in SETS.items():
        for mk,ml in (("full","ARE-full"),("half","ARE-half"),("jaspar","JASPAR")):
            P=C.pairs(tree,mk,"site",keep=S)
            ga=C.group(P,"ra","sub"); gb=C.group(P,"g4","sub")
            if len(ga)<5 or len(gb)<5:
                print(f"  {lab:<16}{ml:<9}{'too few':>19}{'':38}  {len(ga)}/{len(gb)}"); continue
            a=ga[:,0]-ga[:,1]; b=gb[:,0]-gb[:,1]
            ra,pra=p1(a); g4,pg4=p1(b); d,pd=perm(a,b)
            print(f"  {lab:<16}{ml:<9}{ra:>+10.3f} (p{pra:5.3f}){st(pra)}"
                  f"{g4:>+9.3f} (p{pg4:5.3f}){st(pg4)}"
                  f"{d:>+9.3f} (p{pd:5.3f}){st(pd)}  {len(a)}/{len(b)}")
        print()
