#!/usr/bin/env python
"""Trackhub BEDs for AR motif sites in Glut-RA / GABA-4 DARs (lonStrDom2).

One feature per MERGED site, so the two strand reports of one palindromic ARE
collapse into a single feature. Strand and FIMO score come from the
best-scoring hit of that site.

The feature is the BEST-SCORING hit of the site, exactly W bp, so the interval
matches the strand and FIMO score reported with it. phyloP is averaged over
those same W columns, so the displayed interval is the scored interval.

BED9 + 5:  darClass, element, fimoScore, phyloP_oscine, phyloP_wholeTree
score (col 5) is the oscine phyloP rescaled to 0-1000, 500 = zero.

usage: make_ar_track.py <tag>      tag in {full, half, jaspar}
"""
import collections, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
MOT={"full":("motif2",12,"arefull","HOMER de-novo ARE-full"),
     "half":("motif4",8,"arehalf","HOMER de-novo ARE-half"),
     "jaspar":("jaspar",16,"jasparAR","JASPAR MA0007.3 (AR)")}
tag=sys.argv[1]; pref,W,stem,desc=MOT[tag]
RGB={"Glut-RA":"227,26,28","GABA-4":"44,110,127"}
CL={"RA_over_Arco":"Glut-RA","GABA4_over_GABA1":"GABA-4"}
el={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); el[f[3]]=(f[0],int(f[1]),int(f[2]),f[4])
mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
ph={}
for tree in ("subtree222","wholetree"):
    acc=collections.defaultdict(list)
    for l in open(f"{F}/{tree}_ar3_pos.tab"):
        p=l.split("\t")
        if int(p[2])==0: continue
        m,ss,loc,e,site,k=p[0].split("|")
        if m!=tag: continue
        acc[(ss,loc,e,int(site))].append(float(p[5]))
    ph[tree]={k:sum(v)/len(v) for k,v in acc.items() if len(v)>=W/2}
per=collections.defaultdict(list)
for s in ("ra","g4"):
    for i,l in enumerate(open(f"{F}/fimo_{pref}_{s}.txt")):
        if not i: continue
        p=l.rstrip("\n").split("\t"); e=p[2]
        if e not in el: continue
        ch,es,ee,cl=el[e]
        if (cl=="RA_over_Arco")!=(s=="ra"): continue
        st,sp=int(p[3]),int(p[4])
        if sp-st+1!=W: continue
        per[(s,e)].append((es+st-1,es+sp,p[5],float(p[6])))
rows=[]
for (s,e),v in sorted(per.items()):
    ch,_,_,cl=el[e]; loc="MAFB" if e in mafb else "other"; c=CL[cl]
    v.sort(); groups=[]
    for h in v:
        if groups and h[0]<groups[-1][-1][1]: groups[-1].append(h)
        else: groups.append([h])
    for j,g in enumerate(groups):
        a,b,strand,sc=max(g,key=lambda x:x[3])     # the best-scoring hit, exactly W bp
        key=(s,loc,e,j)
        po=ph["subtree222"].get(key); pw=ph["wholetree"].get(key)
        if po is None or pw is None: continue
        score=max(0,min(1000,int(round((po+2.0)/4.0*1000))))
        name=("RA_" if c=="Glut-RA" else "G4_")+f"{e}_{j+1}"
        rows.append((ch,a,b,name,score,strand,a,b,RGB[c],c,e,f"{sc:.3f}",
                     f"{po:.4f}",f"{pw:.4f}",loc))
rows.sort(key=lambda r:(r[0],r[1]))
for sub,fn in ((None,f"tracks/genome_wide/{stem}_sites_lsd2.bed"),("MAFB",f"tracks/mafb_1p5mb/{stem}_sites_MAFB_lsd2.bed")):
    with open(f"{F}/{fn}","w") as o:
        for r in rows:
            if sub and r[14]!=sub: continue
            o.write("\t".join(str(x) for x in r[:14])+"\n")
n=collections.Counter(r[9] for r in rows)
m=collections.Counter(r[9] for r in rows if r[14]=="MAFB")
print(f"  {desc}: {len(rows):,} sites   Glut-RA {n['Glut-RA']:,}   GABA-4 {n['GABA-4']:,}")
print(f"    at MAFB {sum(m.values())}   Glut-RA {m['Glut-RA']}   GABA-4 {m['GABA-4']}")
with open(f"{F}/tracks/{stem}.as","w") as o:
    o.write(f'''table {stem}
"{desc} sites in Glut-RA and GABA-4 DARs (lonStrDom2)"
(
string  chrom;          "Reference sequence chromosome"
uint    chromStart;     "Start position"
uint    chromEnd;       "End position"
string  name;           "Site name"
uint    score;          "Oscine phyloP rescaled 0-1000 (500 = zero)"
char[1] strand;         "Strand of the best-scoring hit"
uint    thickStart;     "Start of thick drawing"
uint    thickEnd;       "End of thick drawing"
uint    reserved;       "Item RGB colour"
string  darClass;       "DAR class: Glut-RA or GABA-4"
string  element;        "Parent DAR"
float   fimoScore;      "FIMO log-odds score of the best hit"
float   phyloPOscine;   "Mean phyloP across the site, oscine subtree"
float   phyloPWholeTree; "Mean phyloP across the site, whole tree"
)
''')
