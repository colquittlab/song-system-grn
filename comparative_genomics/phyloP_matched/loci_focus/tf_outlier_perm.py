#!/usr/bin/env python
"""Westfall-Young max-statistic outlier test on the whole-tree / oscine DoD line.

Null: a motif's cell-type asymmetry is the same as the other regulators' --
i.e. it sits on the trend line. Cell-type labels are permuted within each
motif, carrying each element's oscine AND whole-tree value together, so the
dependence between the two axes is preserved exactly. The whole scatter is
refit on every permutation and the max |studentized residual| recorded, which
controls family-wise error across all motifs without a Bonferroni factor.
"""
import numpy as np, sys, json
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
from scipy import stats
sys.path.insert(0,".")
from tf_core import REG, diffs
NPERM=20000

def load(loc):
    O=diffs("subtree222",loc); W=diffs("wholetree",loc)
    mn=8 if loc else 40; P={}
    for m in REG:
        aO,bO=O.get((m,"RA_over_Arco")),O.get((m,"GABA4_over_GABA1"))
        aW,bW=W.get((m,"RA_over_Arco")),W.get((m,"GABA4_over_GABA1"))
        if aO is None or bO is None or len(aO)<mn or len(bO)<mn: continue
        P[m]=(np.concatenate([aO,bO]),np.concatenate([aW,bW]),len(aO))
    return P

def resids(x,y):
    """leave-one-out studentized residual for every point"""
    n=len(x); t=np.empty(n)
    for i in range(n):
        k=np.ones(n,bool); k[i]=False
        a,b=np.polyfit(x[k],y[k],1)[::-1]
        r=y[k]-(a+b*x[k]); dof=n-3
        sse=(r**2).sum()/dof; mx=x[k].mean(); sxx=((x[k]-mx)**2).sum()
        t[i]=(y[i]-(a+b*x[i]))/np.sqrt(sse*(1+1/(n-1)+(x[i]-mx)**2/sxx))
    return t

def run(loc,label):
    P=load(loc); tfs=list(P); n=len(tfs)
    rng=np.random.default_rng(11)
    x=np.array([P[m][1][:P[m][2]].mean()-P[m][1][P[m][2]:].mean() for m in tfs])
    y=np.array([P[m][0][:P[m][2]].mean()-P[m][0][P[m][2]:].mean() for m in tfs])
    tobs=resids(x,y)
    maxnull=np.empty(NPERM); hits=np.zeros(n)
    for s in range(NPERM):
        xs=np.empty(n); ys=np.empty(n)
        for j,m in enumerate(tfs):
            cO,cW,na=P[m]; idx=rng.permutation(len(cO)); A,B=idx[:na],idx[na:]
            ys[j]=cO[A].mean()-cO[B].mean(); xs[j]=cW[A].mean()-cW[B].mean()
        tn=resids(xs,ys)
        maxnull[s]=np.abs(tn).max()
        hits+=np.abs(tn)>=np.abs(tobs)          # per-motif marginal null too
    pfwer=np.array([(np.sum(maxnull>=abs(t))+1)/(NPERM+1) for t in tobs])
    pmarg=(hits+1)/(NPERM+1)
    print("\n=== %s   n=%d  (%d permutations)"%(label,n,NPERM))
    print("    TF      reg      whole    oscine   stud.resid |  p marginal   p FWER")
    o=np.argsort(pfwer)
    for i in o:
        print("    %-7s %-8s %+.3f  %+.3f  %+6.2f     |  %.4f       %.4f %s"%(
            tfs[i],REG[tfs[i]],x[i],y[i],tobs[i],pmarg[i],pfwer[i],
            "*" if pfwer[i]<.05 else ""))
    json.dump([dict(tf=tfs[i],reg=REG[tfs[i]],x=float(x[i]),y=float(y[i]),
                    t=float(tobs[i]),p_marg=float(pmarg[i]),p_fwer=float(pfwer[i]))
               for i in o],open(result(f"tf_outlier_perm_{loc or 'gw'}.json", loc),"w"),indent=1)

run(None,"genome-wide")
run("MAFB","at MAFB (gene +/-1.5 Mb)")
run("domain","in the DHX35-TOP1 domain")
