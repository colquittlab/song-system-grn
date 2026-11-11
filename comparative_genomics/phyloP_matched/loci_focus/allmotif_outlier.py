#!/usr/bin/env python
"""Are the eRegulon TFs outliers from the whole-tree / oscine DoD trend?

The 182 motif archetypes supply the expectation: they fix the trend line and,
more importantly, the spread of residuals around it. Each eRegulon TF is then
tested against that background. Three readings, from least to most assumption-
laden:
  empirical  -- fraction of archetypes further from the line than this TF
  robust z   -- residual / MAD-scaled background spread
  meta z     -- residual / sqrt(own sampling var + between-motif tau^2)
The last is the only one that knows some motifs are measured far better than
others. AR is measured with the SCENIC+ matrix for its cistrome (from
metacluster_185.1), like every other regulator; the metacluster_185.1 archetype
is left out of the fit when AR is the point being tested. AR* is the HOMER
ARE-full motif per instance, a labelled reference.
"""
import numpy as np, json, sys
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
from scipy import stats
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0,F)
from tf_core import REG, SAME_AS
AR_ARCH="metacluster_185.1.jaspar__MA0007.3"
NAME={"KLF12":"KLF"}; DROP={"KLF15"}|set(SAME_AS)

def wls(x,y,w):
    W=w.sum(); mx=(w*x).sum()/W; my=(w*y).sum()/W
    b=(w*(x-mx)*(y-my)).sum()/(w*(x-mx)**2).sum(); return my-b*mx,b

def paule_mandel(x,y,s2):
    n=len(x); lo,hi=0.0,max(1e-12,np.var(y,ddof=1)*50)
    def Q(t2):
        w=1/(s2+t2); a,b=wls(x,y,w); r=y-(a+b*x)
        return (r**2/(s2+t2)).sum()-(n-2)
    if Q(lo)<=0: return 0.0
    for _ in range(200):
        m=(lo+hi)/2
        if Q(m)>0: lo=m
        else: hi=m
    return (lo+hi)/2

def bh(ps):
    ps=np.array(ps); o=np.argsort(ps); q=ps[o]*len(ps)/np.arange(1,len(ps)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty(len(ps)); r[o]=np.clip(q,0,1); return r

def run(tagl,mn,label,suffix=""):
    rows=json.load(open(f"{F}/results/compare_mafb_spans/pairs_all.json"))
    arch={r["motif"]:r[tagl] for r in rows if r["set"]=="arch" and r.get(tagl)}
    ereg={r["motif"]:r[tagl] for r in rows if r["set"]=="ereg" and r.get(tagl)}
    arch={k:v for k,v in arch.items() if v["na"]>=mn and v["nb"]>=mn}
    ereg={NAME.get(k,k):v for k,v in ereg.items()
          if k not in DROP and v["na"]>=mn and v["nb"]>=mn}
    # the manuscript's AR PWM, scored by paired_core -- a different instance set
    from tf_core import diffs as _d
    try:
        O=_d("subtree222",None if tagl=="gw" else tagl)
        W=_d("wholetree", None if tagl=="gw" else tagl)
        aO,bO=O[("AR*","RA_over_Arco")],O[("AR*","GABA4_over_GABA1")]
        aW,bW=W[("AR*","RA_over_Arco")],W[("AR*","GABA4_over_GABA1")]
        assert len(aO)==len(aW) and len(bO)==len(bW)
        na,nb=len(aO),len(bO)
        if na>=mn and nb>=mn:
            ereg["AR*"]=dict(na=na,nb=nb,
                y=float(aO.mean()-bO.mean()), x=float(aW.mean()-bW.mean()),
                vy=float(aO.var(ddof=1)/na+bO.var(ddof=1)/nb),
                vx=float(aW.var(ddof=1)/na+bW.var(ddof=1)/nb),
                cxy=float(np.cov(aO,aW)[0,1]/na+np.cov(bO,bW)[0,1]/nb))
    except Exception as e:
        print("   (AR* unavailable: %s)"%e)
    names=list(arch); X=np.array([arch[m]["x"] for m in names]); Y=np.array([arch[m]["y"] for m in names])
    b,a=np.polyfit(X,Y,1); r_,p_=stats.pearsonr(X,Y)
    S2=np.array([arch[m]["vy"]+b*b*arch[m]["vx"]-2*b*arch[m]["cxy"] for m in names])
    tau2=paule_mandel(X,Y,S2)
    R=Y-(a+b*X); sig=1.4826*np.median(np.abs(R-np.median(R)))
    print(f"\n{'='*96}\n=== {label}:  {len(names)} archetypes as background, {len(ereg)} eRegulon TFs tested")
    print(f"    trend  y = {a:+.4f} {b:+.4f} x     r = {r_:+.3f} (p = {p_:.2g})")
    print(f"    residual spread: robust sd {sig:.4f}   median sampling SE {np.sqrt(np.median(S2)):.4f}"
          f"   tau {np.sqrt(tau2):.4f}")
    frac=tau2/(tau2+np.median(S2))
    print(f"    -> {100*frac:.0f}% of the scatter around the line is real between-motif variation")
    out=[]
    for m,d in ereg.items():
        keep=[i for i,nm in enumerate(names) if nm!=(AR_ARCH if m=="AR" else None)]
        xk,yk,s2k=X[keep],Y[keep],S2[keep]
        bb,aa=np.polyfit(xk,yk,1)
        rk=yk-(aa+bb*xk); sg=1.4826*np.median(np.abs(rk-np.median(rk)))
        res=d["y"]-(aa+bb*d["x"])
        zr=res/sg; p_emp=(np.sum(np.abs(rk)>=abs(res))+1)/(len(rk)+1)
        t2=paule_mandel(xk,yk,s2k)
        w=1/(s2k+t2); aw,bw=wls(xk,yk,w)
        W=w.sum(); mxw=(w*xk).sum()/W
        vline=1/W+(d["x"]-mxw)**2/(w*(xk-mxw)**2).sum()
        s2i=d["vy"]+bw*bw*d["vx"]-2*bw*d["cxy"]
        zm=(d["y"]-(aw+bw*d["x"]))/np.sqrt(s2i+t2+vline)
        out.append(dict(tf=m,reg=REG.get(m.rstrip("*"),"?"),x=d["x"],y=d["y"],res=res,se=np.sqrt(s2i),
                        z_rob=zr,p_emp=p_emp,z_meta=zm,p_meta=2*stats.norm.sf(abs(zm)),
                        na=d["na"],nb=d["nb"]))
    for k,lab in (("p_emp","q_emp"),("p_meta","q_meta")):
        for r,q in zip(out,bh([x[k] for x in out])): r[lab]=q
    out.sort(key=lambda r:r["p_meta"])
    print(f"\n    {'TF':<7}{'regulon':<9}{'whole':>8}{'oscine':>8}{'resid':>8}{'ownSE':>8}"
          f"{'z_rob':>7}{'p_emp':>8}{'q_emp':>7}{'z_meta':>8}{'p_meta':>8}{'q_meta':>8}  n")
    for r in out:
        st="*" if r["q_meta"]<.05 else (":" if r["p_meta"]<.05 else " ")
        print(f"    {r['tf']:<7}{r['reg']:<9}{r['x']:>+8.3f}{r['y']:>+8.3f}{r['res']:>+8.3f}"
              f"{r['se']:>8.4f}{r['z_rob']:>+7.2f}{r['p_emp']:>8.4f}{r['q_emp']:>7.3f}"
              f"{r['z_meta']:>+8.2f}{r['p_meta']:>8.4f}{r['q_meta']:>8.3f}{st} {r['na']}/{r['nb']}")
    json.dump(dict(fit=dict(a=a,b=b,r=r_,p=p_,sig=sig,tau=float(np.sqrt(tau2))),
                   arch=[dict(motif=m,**arch[m]) for m in names],ereg=out),
              open(result(f"allmotif_outlier_{tagl}{suffix}.json", tagl),"w"),indent=1,default=float)
    return out

run("gw",40,"genome-wide")
try: run("MAFB",8,"at the MAFB locus (gene +/-1.5 Mb)")
except Exception as e: print("\n  MAFB: %s"%e)
# the domain is small (118 DARs): >=8 per class matches the +/-1.5 Mb run; >=4 is a
# sensitivity run that admits the AR archetype (4 Glut-RA / 6 GABA-4 elements)
run("domain",8,"in the DHX35-TOP1 domain, >=8 elements per class")
run("domain",4,"in the DHX35-TOP1 domain, >=4 elements per class (sensitivity)","_min4")
