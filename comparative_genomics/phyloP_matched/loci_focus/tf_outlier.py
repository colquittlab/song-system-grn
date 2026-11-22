#!/usr/bin/env python
"""Outlier tests on the whole-tree / oscine DoD trend line.

Two tests per locus set:
  (1) classic externally-studentized (deleted) residual, Bonferroni over all motifs
  (2) meta-regression residual, using each motif's own sampling SE plus an
      estimated between-motif dispersion tau
Both refit the line with the tested point removed, so a point cannot hide
inside a line it defined.
"""
import numpy as np, sys, json
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
from scipy import stats
sys.path.insert(0,".")
from tf_core import REG, diffs

def moments(loc):
    """per-motif DoD on each tree with sampling var/cov, from the same elements"""
    O=diffs("subtree222",loc); W=diffs("wholetree",loc)
    mn=8 if loc else 40; out={}
    for m in REG:
        aO=O.get((m,"RA_over_Arco")); bO=O.get((m,"GABA4_over_GABA1"))
        aW=W.get((m,"RA_over_Arco")); bW=W.get((m,"GABA4_over_GABA1"))
        if aO is None or bO is None or aW is None or bW is None: continue
        if len(aO)<mn or len(bO)<mn: continue
        assert len(aO)==len(aW) and len(bO)==len(bW), m
        na,nb=len(aO),len(bO)
        y=aO.mean()-bO.mean(); x=aW.mean()-bW.mean()
        vy=aO.var(ddof=1)/na + bO.var(ddof=1)/nb
        vx=aW.var(ddof=1)/na + bW.var(ddof=1)/nb
        cxy=np.cov(aO,aW)[0,1]/na + np.cov(bO,bW)[0,1]/nb
        out[m]=dict(x=x,y=y,vx=vx,vy=vy,cxy=cxy,na=na,nb=nb)
    return out

def wls(x,y,w):
    W=w.sum(); mx=(w*x).sum()/W; my=(w*y).sum()/W
    b=(w*(x-mx)*(y-my)).sum()/(w*(x-mx)**2).sum()
    return my-b*mx, b

def paule_mandel(x,y,s2,lo=0.0,hi=None):
    """tau^2 solving sum r^2/(s^2+tau^2) = n-2, with WLS refit each step"""
    n=len(x)
    if hi is None: hi=max(1e-12,np.var(y,ddof=1)*50)
    def Q(t2):
        w=1/(s2+t2); a,b=wls(x,y,w); r=y-(a+b*x)
        return (r**2/(s2+t2)).sum()-(n-2)
    if Q(lo)<=0: return 0.0
    for _ in range(200):
        mid=(lo+hi)/2
        if Q(mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2

def run(loc,label):
    D=moments(loc); tfs=list(D); n=len(tfs)
    x=np.array([D[t]["x"] for t in tfs]); y=np.array([D[t]["y"] for t in tfs])
    a,b=np.polyfit(x,y,1)[::-1]
    r,p=stats.pearsonr(x,y)
    # sampling variance of the vertical residual, given slope b
    s2=np.array([D[t]["vy"]+b*b*D[t]["vx"]-2*b*D[t]["cxy"] for t in tfs])
    tau2=paule_mandel(x,y,s2)
    print("\n=== %s   n=%d   slope %+.4f   r %+.3f (p %.3f)"%(label,n,b,r,p))
    print("    median sampling SE %.4f   between-motif tau %.4f  (ratio %.2f)"%(
        np.sqrt(np.median(s2)),np.sqrt(tau2),np.sqrt(tau2/np.median(s2))))
    rows=[]
    for i,t in enumerate(tfs):
        k=np.ones(n,bool); k[i]=False
        # (1) classic deleted studentized residual, OLS on the other n-1
        ai,bi=np.polyfit(x[k],y[k],1)[::-1]
        res=y[k]-(ai+bi*x[k]); dof=n-1-2
        sse=(res**2).sum()/dof
        mx=x[k].mean(); sxx=((x[k]-mx)**2).sum()
        se_pred=np.sqrt(sse*(1+1/(n-1)+(x[i]-mx)**2/sxx))
        tstat=(y[i]-(ai+bi*x[i]))/se_pred
        p_ols=2*stats.t.sf(abs(tstat),dof)
        # (2) meta-regression: weighted refit without i, residual vs own SE + tau
        s2k=np.array([D[u]["vy"]+bi*bi*D[u]["vx"]-2*bi*D[u]["cxy"] for u in tfs])
        t2=paule_mandel(x[k],y[k],s2k[k])
        w=1/(s2k[k]+t2); aw,bw=wls(x[k],y[k],w)
        # predictive SE of the line at x_i
        W=w.sum(); mxw=(w*x[k]).sum()/W
        var_line=1/W + (x[i]-mxw)**2/(w*(x[k]-mxw)**2).sum()
        s2i=D[t]["vy"]+bw*bw*D[t]["vx"]-2*bw*D[t]["cxy"]
        z=(y[i]-(aw+bw*x[i]))/np.sqrt(s2i+t2+var_line)
        p_meta=2*stats.norm.sf(abs(z))
        rows.append(dict(tf=t,reg=REG[t],x=x[i],y=y[i],resid=y[i]-(ai+bi*x[i]),
                         se=np.sqrt(s2i),t=tstat,p_ols=p_ols,z=z,p_meta=p_meta))
    rows.sort(key=lambda d:d["p_meta"])
    print("    TF      reg      whole    oscine   resid    ownSE  |  t      p_ols  bonf  |  z      p_meta bonf")
    for d in rows:
        print("    %-7s %-8s %+.3f  %+.3f  %+.3f  %.4f | %+5.2f  %.4f %.3f | %+5.2f  %.4f %.3f"%(
            d["tf"],d["reg"],d["x"],d["y"],d["resid"],d["se"],
            d["t"],d["p_ols"],min(1,d["p_ols"]*n),d["z"],d["p_meta"],min(1,d["p_meta"]*n)))
    json.dump(rows,open(result(f"tf_outlier_{loc or 'gw'}.json", loc),"w"),indent=1,default=float)
    return rows

run(None,"genome-wide")
run("MAFB","at MAFB (gene +/-1.5 Mb)")
run("domain","in the DHX35-TOP1 domain")
