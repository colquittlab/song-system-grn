#!/usr/bin/env python
"""How much does the outgroup -> oscine change alter the PWM match?

For each AR site, build the outgroup consensus and the oscine consensus from
the alignment column majorities, score both as PWM log-odds against a uniform
background, and take the difference. Positive = the oscine sequence is a
BETTER match to the motif than the outgroup sequence; negative = worse.

This replaces the binary switch call with a graded, functional quantity: a
change at a high-information column moves the score a lot, one at a degenerate
column barely at all.
"""
import json, sys, collections, numpy as np
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
D=f"{F}/maf_switch"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC=0.001                                   # pseudocount, keeps log finite
def homer_pwm(f): return np.array([[float(x) for x in l.split()] for l in open(f)
                                   if not l.startswith(">")])
def meme_pwm(f,W):
    b=open(f).read().split("letter-probability")[1].split("\n")[1:]
    return np.array([[float(x) for x in l.split()] for l in b if len(l.split())==4][:W])
PWM={"full":homer_pwm(f"{HM}/motif2.motif"),"jaspar":meme_pwm(f"{F}/ar_jaspar.meme",16)}
LO={k:np.log2((v+PC)/(1+4*PC)/0.25) for k,v in PWM.items()}   # log-odds vs uniform
IC={k:np.array([2+(p*np.log2(np.clip(p,1e-9,None))).sum() for p in v]) for k,v in PWM.items()}
B={"A":0,"C":1,"G":2,"T":3}
DOM={l.strip() for l in open(f"{D}/data/domain_dars.txt") if l.strip()}
parent=lambda s:s.split("_",1)[1].rsplit("_",1)[0]
cls=lambda s:"Glut-RA" if s.startswith("RA_") else "GABA-4"

def deltas(path,minO,minG):
    """-> list of dicts, one per site that is fully callable in both groups"""
    out=[]
    for r in json.load(open(path)):
        mo=r["motif"]
        if mo not in LO: continue
        W=len(LO[mo])
        if len(r["cols"])!=W: continue
        so=sg=0.0; ok=True; ch=[]
        for c in r["cols"]:
            if c["oscN"]<minO or c["outN"]<minG or c["osc"] not in B or c["out"] not in B:
                ok=False; break
            k=c["k"]-1
            so+=LO[mo][k][B[c["osc"]]]; sg+=LO[mo][k][B[c["out"]]]
            if c["osc"]!=c["out"]:
                ch.append((c["k"],c["out"],c["osc"],
                           LO[mo][k][B[c["osc"]]]-LO[mo][k][B[c["out"]]],IC[mo][k]))
        if not ok: continue
        out.append(dict(motif=mo,site=r["site"],cls=cls(r["site"]),
                        dom=parent(r["site"]) in DOM,osc=so,out=sg,delta=so-sg,changes=ch))
    return out

rng=np.random.default_rng(0)
def perm(a,b,n=20000):
    o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0
    for _ in range(n):
        p=rng.permutation(c); h+=abs(p[:k].mean()-p[k:].mean())>=abs(o)
    return o,(h+1)/(n+1)
def p1(a,n=20000):
    s=rng.choice([-1,1],(n,len(a)))
    return a.mean(),((np.abs((s*a).mean(1))>=abs(a.mean())).sum()+1)/(n+1)

if __name__=="__main__":
    which=sys.argv[1] if len(sys.argv)>1 else "363"
    path,mo,mg,lab={"363":(f"{D}/data/sites363.json",40,60,"avian363 (146 oscine / 217 outgroup)"),
                    "74":(f"{D}/data74/sites74.json",6,15,"74-way (17 / 56)")}[which]
    R=deltas(path,mo,mg)
    print(f"=== {lab}: PWM log-odds, oscine consensus minus outgroup consensus\n")
    st=lambda p:"*" if p<.05 else " "
    print(f"  {'scope':<22}{'motif':<9}{'Glut-RA mean':>22}{'GABA-4 mean':>22}{'difference':>20}  n")
    for scope,sel in (("MAFB +/-1.5 Mb",lambda x:True),("DHX35-TOP1 domain",lambda x:x["dom"])):
        for mk,ml in (("full","ARE-full"),("jaspar","JASPAR")):
            a=np.array([x["delta"] for x in R if x["motif"]==mk and sel(x) and x["cls"]=="Glut-RA"])
            b=np.array([x["delta"] for x in R if x["motif"]==mk and sel(x) and x["cls"]=="GABA-4"])
            if len(a)<3 or len(b)<3:
                print(f"  {scope:<22}{ml:<9}{'too few':>22}  n={len(a)}/{len(b)}"); continue
            ma,pa=p1(a); mb,pb=p1(b); d,pd=perm(a,b)
            print(f"  {scope:<22}{ml:<9}{ma:>+13.3f} (p{pa:5.3f}){st(pa)}"
                  f"{mb:>+13.3f} (p{pb:5.3f}){st(pb)}{d:>+11.3f} (p{pd:5.3f}){st(pd)}"
                  f"  {len(a)}/{len(b)}")
        print()
    json.dump(R,open(f"{D}/data/pwm_delta_{which}.json","w"),indent=1,default=float)
    print(f"  wrote data/pwm_delta_{which}.json  ({len(R)} scorable sites)")
