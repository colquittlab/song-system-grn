#!/usr/bin/env python
"""Trackhub BED for the HOMER de-novo ARE-full sites in Glut-RA / GABA-4 DARs.

One feature per MERGED site (overlapping FIMO hits -- the two strand reports of
one ARE -- collapse into a single feature). Strand and FIMO score come from the
best-scoring hit of that site. Carries per-site phyloP on both trees so the
track can be shaded and inspected directly in the browser.

BED9 + 5:  class, element, fimoScore, phyloP_oscine, phyloP_wholetree
score (col 5) is the oscine phyloP rescaled to 0-1000, 500 = zero.
"""
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
W=12
RGB={"Glut-RA":"227,26,28","GABA-4":"44,110,127"}
el={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); el[f[3]]=(f[0],int(f[1]),int(f[2]),f[4])
mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
CL={"RA_over_Arco":"Glut-RA","GABA4_over_GABA1":"GABA-4"}
# phyloP per site, keyed motif|class|locus|element|index
ph={}
for tree in ("subtree222","wholetree"):
    d={}
    for l in open(f"{F}/{tree}_ar3_sites.tab"):
        p=l.split("\t")
        if int(p[2])==0: continue
        d[p[0]]=float(p[5])
    ph[tree]=d
rows=[]
per=collections.defaultdict(list)
for s in ("ra","g4"):
    for i,l in enumerate(open(f"{F}/fimo_motif2_{s}.txt")):
        if not i: continue
        p=l.rstrip("\n").split("\t"); e=p[2]
        if e not in el: continue
        ch,es,ee,cl=el[e]
        if (cl=="RA_over_Arco")!=(s=="ra"): continue
        st,sp=int(p[3]),int(p[4])
        if sp-st+1!=W: continue
        per[(s,e)].append((es+st-1,es+sp,p[5],float(p[6])))
for (s,e),v in sorted(per.items()):
    ch,_,_,cl=el[e]; loc="MAFB" if e in mafb else "other"; c=CL[cl]
    v.sort(); groups=[]
    for h in v:
        if groups and h[0]<groups[-1][-1][1]: groups[-1].append(h)
        else: groups.append([h])
    for j,g in enumerate(groups):
        a=min(x[0] for x in g); b=max(x[1] for x in g)
        _,_,strand,sc=max(g,key=lambda x:x[3])
        key=f"full|{s}|{loc}|{e}|{j}"
        po=ph["subtree222"].get(key); pw=ph["wholetree"].get(key)
        if po is None or pw is None: continue
        score=max(0,min(1000,int(round((po+2.0)/4.0*1000))))
        name=("RA_" if c=="Glut-RA" else "G4_")+f"{e}_{j+1}"
        rows.append((ch,a,b,name,score,strand,a,b,RGB[c],c,e,f"{sc:.3f}",
                     f"{po:.4f}",f"{pw:.4f}",loc))
rows.sort(key=lambda r:(r[0],r[1]))
with open(f"{F}/tracks/genome_wide/arefull_sites_lsd2.bed","w") as o:
    for r in rows: o.write("\t".join(str(x) for x in r[:14])+"\n")
MAFBN={r[3] for r in rows if r[14]=="MAFB"}
with open(f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","w") as o:
    for r in rows:
        if r[14]=="MAFB": o.write("\t".join(str(x) for x in r[:14])+"\n")
n=collections.Counter(r[9] for r in rows)
m=collections.Counter(r[9] for r in rows if r[14]=="MAFB")
print(f"  total sites {len(rows):,}   Glut-RA {n['Glut-RA']:,}   GABA-4 {n['GABA-4']:,}")
print(f"  at MAFB      {sum(m.values()):,}   Glut-RA {m['Glut-RA']}   GABA-4 {m['GABA-4']}")
