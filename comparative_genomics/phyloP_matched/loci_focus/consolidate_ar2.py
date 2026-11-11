#!/usr/bin/env python
"""Parts 2 and 3 of the AR report: switch + motif scoring, and the cross-motif set."""
import json, collections, sys, numpy as np
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
sys.path.insert(0,"."); sys.path.insert(0,"maf_switch")
from ar_highic import IC
M="maf_switch"
rng=np.random.default_rng(0); NP=20000
DOM={l.strip() for l in open(f"{M}/data/domain_dars.txt") if l.strip()}
parent=lambda s:s.split("_",1)[1].rsplit("_",1)[0]
OUT={}
# --- 2a: switch rates + permutation, by alignment x scope x motif
def rows(path,mo,mg,restrict):
    out=collections.defaultdict(list)
    for r in json.load(open(path)):
        if restrict and parent(r["site"]) not in DOM: continue
        c="Glut-RA" if r["site"].startswith("RA_") else "GABA-4"
        te=sw=0
        for x in r["cols"]:
            if x["oscN"]<mo or x["outN"]<mg: continue
            te+=1
            if x["osc"]!=x["out"] and x["oscF"]>=.8 and x["outF"]>=.8: sw+=1
        out[r["motif"]].append((c,sw,te))
    return out
def ptest(rs):
    cl=np.array([x[0]=="Glut-RA" for x in rs])
    sw=np.array([x[1] for x in rs],float); te=np.array([x[2] for x in rs],float)
    f=lambda mk:(sw[mk].sum()/te[mk].sum()-sw[~mk].sum()/te[~mk].sum())
    obs=f(cl); n=len(rs); k=int(cl.sum()); hit=0; r2=np.random.default_rng(0)
    for _ in range(NP):
        p=np.zeros(n,bool); p[r2.choice(n,k,replace=False)]=True
        if abs(f(p))>=abs(obs)-1e-12: hit+=1
    return dict(diff=float(obs),p=float((hit+1)/(NP+1)),
                sw_ra=float(sw[cl].sum()),te_ra=float(te[cl].sum()),
                sw_g4=float(sw[~cl].sum()),te_g4=float(te[~cl].sum()),
                sites_ra=int(cl.sum()),sites_g4=int((~cl).sum()))
OUT["switch"]={}
for al,path,mo,mg,sp in (("avian363",f"{M}/data/sites363.json",40,60,"146 / 217"),
                         ("74way",f"{M}/data74/sites74.json",6,15,"17 / 56")):
    OUT.setdefault("split",{})[al]=sp
    for scope,res in (("1.5Mb",False),("domain",True)):
        D=rows(path,mo,mg,res)
        for mk in ("full","jaspar"):
            if len(D.get(mk,[]))>=4: OUT["switch"][f"{al}|{scope}|{mk}"]=ptest(D[mk])
# --- 2b: IC and phyloP at switch vs non-switch columns (avian363)
PH={}
for tree in ("subtree222","wholetree"):
    dd={}
    for l in open(f"{tree}_ar3_pos.tab"):
        p=l.split("\t"); mm,s,loc,e,site,k=p[0].split("|")
        if int(p[2])==0: continue
        dd[(mm,e,site,int(k))]=float(p[5])
    PH[tree]=dd
idx={}
for bed,mo in (("tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","full"),("tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed","jaspar")):
    per=collections.defaultdict(int)
    for l in open(bed):
        x=l.split("\t"); e=x[10]; idx[(mo,x[3])]=(e,str(per[e])); per[e]+=1
SW=collections.defaultdict(list); NS=collections.defaultdict(list)
ICs=collections.defaultdict(list); ICn=collections.defaultdict(list)
for r in json.load(open(f"{M}/data/sites363.json")):
    key=idx.get((r["motif"],r["site"]))
    if key is None: continue
    e,si=key
    for c in r["cols"]:
        if c["oscN"]<40 or c["outN"]<60: continue
        k=c["k"]-1
        sw=c["osc"]!=c["out"] and c["oscF"]>=.8 and c["outF"]>=.8
        (ICs if sw else ICn)[r["motif"]].append(float(IC[r["motif"]][k]))
        for tree in ("subtree222","wholetree"):
            v=PH[tree].get((r["motif"],e,si,k))
            if v is not None: (SW if sw else NS)[(r["motif"],tree)].append(v)
OUT["switch_phylop"]={}
for mo in ("full","jaspar"):
    for tree in ("subtree222","wholetree"):
        a=np.array(SW[(mo,tree)]); b=np.array(NS[(mo,tree)])
        if len(a)<3: continue
        bs=np.array([np.mean(rng.choice(a,len(a)))-np.mean(rng.choice(b,len(b))) for _ in range(NP)])
        p=max(2*min(float((bs<=0).mean()),float((bs>=0).mean())),1/NP)
        OUT["switch_phylop"][f"{mo}|{tree}"]=dict(sw=float(a.mean()),n_sw=len(a),
                                                  ns=float(b.mean()),n_ns=len(b),
                                                  diff=float(a.mean()-b.mean()),p=p)
OUT["switch_ic"]={mo:dict(sw=float(np.mean(ICs[mo])),n_sw=len(ICs[mo]),
                          ns=float(np.mean(ICn[mo])),n_ns=len(ICn[mo]),
                          hi_sw=float(np.mean([x>=1 for x in ICs[mo]])),
                          hi_ns=float(np.mean([x>=1 for x in ICn[mo]])))
                 for mo in ("full","jaspar")}
# --- 2c: PWM match change, with the within-element control
R=json.load(open(f"{M}/data/pwm_control.json"))
def p1(a):
    s=rng.choice([-1,1],(NP,len(a)))
    return float(a.mean()),float(((np.abs((s*a).mean(1))>=abs(a.mean())).sum()+1)/(NP+1))
def perm(a,b):
    o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0; r2=np.random.default_rng(0)
    for _ in range(NP):
        p=r2.permutation(c); h+=abs(p[:k].mean()-p[k:].mean())>=abs(o)
    return float(o),float((h+1)/(NP+1))
OUT["pwm"]={}
for scope,sel in (("1.5Mb",lambda r:True),("domain",lambda r:r["dom"])):
    for mk in ("full","jaspar"):
        cellv={}
        for c in ("Glut-RA","GABA-4"):
            S=[r for r in R if r["motif"]==mk and r["cls"]==c and sel(r)]
            if len(S)<3: continue
            a=np.array([r["paired"] for r in S]); mm,pp=p1(a)
            cellv[c]=dict(site=float(np.mean([r["site_delta"] for r in S])),
                          ctrl=float(np.mean([r["ctrl_mean"] for r in S])),
                          paired=mm,p=pp,n=len(S))
        if len(cellv)==2:
            A=np.array([r["paired"] for r in R if r["motif"]==mk and r["cls"]=="Glut-RA" and sel(r)])
            B=np.array([r["paired"] for r in R if r["motif"]==mk and r["cls"]=="GABA-4" and sel(r)])
            dd,pd=perm(A,B); cellv["dod"]=dict(d=dd,p=pd)
            OUT["pwm"][f"{scope}|{mk}"]=cellv
# --- 3: the SCENIC+ MAFB regulator set, and the archetype background
OUT["tf"]={k:json.load(open(result(f"tf_tab_{k}.json", k.rsplit("_", 1)[1]))) for k in
           ("subtree222_gw","wholetree_gw","subtree222_MAFB","subtree222_domain")}
OUT["outlier"]={k:json.load(open(result(f"allmotif_outlier_{k}.json", k))) for k in ("gw","MAFB","domain")}
json.dump(OUT,open("results/compare_mafb_spans/ar_report_part23.json","w"),indent=1)
print("  switch cells:",len(OUT["switch"]),"  pwm cells:",len(OUT["pwm"]),
      "  TF tables:",{k:len(v) for k,v in OUT["tf"].items()})
