#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 48G
#SBATCH -J astro
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/astro-%j.out
set -eu
export PATH=/private/groups/colquittlab/opt/ucsc_utilities:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; cd $F
cat e_ra.bed e_ast.bed | cut -f1-4 | sort -k1,1 -k2,2n > as_dars.bed
echo "  DARs: $(wc -l < as_dars.bed)"
python3 - <<'PY'
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
coord={}
for s in ("ra","ast"):
    for l in open(f"{F}/e_{s}.bed"):
        f=l.rstrip("\n").split("\t"); coord[f[3]]=(f[0],int(f[1]))
iv=collections.defaultdict(list); n=0
for tf in ("Fosl2","Mafb","Nfia"):
    for s in ("ra","ast"):
        for i,l in enumerate(open(f"{F}/efimo_{tf}_{s}.txt")):
            if not i or not l.strip(): continue
            p=l.rstrip("\n").split("\t")
            if len(p)<6 or p[2] not in coord: continue
            iv[(tf,p[2])].append((int(p[3])-1,int(p[4]))); n+=1
o=open(f"{F}/as_merged.bed","w"); k=0
for (m,e),v in iv.items():
    v.sort(); cur=list(v[0]); mg=[]
    for a,b in v[1:]:
        if a<=cur[1]: cur[1]=max(cur[1],b)
        else: mg.append(tuple(cur)); cur=[a,b]
    mg.append(tuple(cur)); c,s0=coord[e]
    for a,b in mg: k+=1; o.write(f"{c}\t{s0+a}\t{s0+b}\t{m}@{e}@{k}\n")
o.close(); print(f"  {n} hits -> {k} merged over {len(iv)} pairs")
PY
sort -k1,1 -k2,2n as_merged.bed -o as_merged.bed
for t in subtree222 wholetree; do
  bigWigAverageOverBed $H/${t}_perbase_lsd2.bw as_dars.bed   ${t}_as_dars.tab
  bigWigAverageOverBed $H/${t}_perbase_lsd2.bw as_merged.bed ${t}_as_hits.tab
done
echo ASTRODONE
