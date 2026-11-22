import numpy as np, sys, json
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
sys.path.insert(0,".")
from tf_core import REG, diffs, bh
rng=np.random.default_rng(1)
def perm(a,b,n=20000):
    o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0
    for _ in range(n):
        p=rng.permutation(c); h+= abs(p[:k].mean()-p[k:].mean())>=abs(o)
    return o,(h+1)/(n+1)
def p1(a,n=20000):   # one-sample sign-flip
    o=a.mean(); s=rng.choice([-1,1],(n,len(a)))
    return o,((np.abs((s*a).mean(1))>=abs(o)).sum()+1)/(n+1)
for tree in ("subtree222","wholetree"):
  for loc,mn in ((None,40),("MAFB",8),("domain",8)):
    D=diffs(tree,loc); rows=[]
    for m in REG:
        a=D.get((m,"RA_over_Arco"),np.array([])); b=D.get((m,"GABA4_over_GABA1"),np.array([]))
        if len(a)<mn or len(b)<mn: continue
        ra,pra=p1(a); g4,pg4=p1(b); d,pd=perm(a,b)
        rows.append(dict(tf=m,reg=REG[m],ra=ra,pra=pra,g4=g4,pg4=pg4,dod=d,p=pd,na=len(a),nb=len(b)))
    for k in ("pra","pg4","p"):
        for r,q in zip(rows,bh([r[k] for r in rows])): r["q_"+k]=q
    rows.sort(key=lambda r:-r["ra"])
    print(f"\n=== {tree} {loc or 'genome-wide'}  ({len(rows)} TFs)")
    print("  TF      reg      Glut-RA   p      q    | GABA-4    p      q    |  DoD      p      q    |  n")
    for r in rows:
        print("  %-7s %-8s %+.3f %6.4f %6.4f | %+.3f %6.4f %6.4f | %+.3f %6.4f %6.4f | %d/%d"%(
            r["tf"],r["reg"],r["ra"],r["pra"],r["q_pra"],r["g4"],r["pg4"],r["q_pg4"],
            r["dod"],r["p"],r["q_p"],r["na"],r["nb"]))
    json.dump(rows,open(result(f"tf_tab_{tree}_{loc or 'gw'}.json", loc),"w"),indent=1)
