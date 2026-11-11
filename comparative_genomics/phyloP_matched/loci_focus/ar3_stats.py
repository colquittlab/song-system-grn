import sys,json,numpy as np; sys.path.insert(0,".")
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
st=lambda p:"*" if p<.05 else " "
res=[]
for tree,tl in (("subtree222","oscine"),("wholetree","whole")):
    print("\n=== %s tree: motif minus rest of its own enhancer"%tl)
    print("  %-8s %-9s %-11s %19s %19s %19s   n"%("motif","unit","scope","Glut-RA","GABA-4","DoD"))
    for m in ("full","half","jaspar"):
        for u in ("instance","site","element"):
            P=C.pairs(tree,m,u)
            for loc,lab in ((None,"genome-wide"),("MAFB","at MAFB")):
                ga=C.group(P,"ra",loc); gb=C.group(P,"g4",loc)
                if len(ga)<5 or len(gb)<5: continue
                a=ga[:,0]-ga[:,1]; b=gb[:,0]-gb[:,1]
                ra,pra=p1(a); g4,pg4=p1(b); d,pd=perm(a,b)
                print("  %-8s %-9s %-11s  %+.3f (p %.4f)%s %+.3f (p %.4f)%s %+.3f (p %.4f)%s %d/%d"%(
                    m,u,lab,ra,pra,st(pra),g4,pg4,st(pg4),d,pd,st(pd),len(a),len(b)))
                res.append(dict(tree=tl,motif=m,unit=u,scope=lab,ra=ra,pra=pra,
                                g4=g4,pg4=pg4,dod=d,p=pd,na=len(a),nb=len(b)))
        print()
json.dump(res,open("results/compare_gw_mafb1p5mb/ar3_stats.json","w"),indent=1,default=float)
