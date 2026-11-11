#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 48G
#SBATCH -c 8
#SBATCH -J tfsc
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/tfsc-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/ucsc_utilities:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; cd $F
fimo --verbosity 1 --thresh 1e-4 --text tf_scenic.meme am_dars.fa > tfsc_hits.txt
echo "  hits: $(( $(wc -l < tfsc_hits.txt) - 1 ))"
python3 - <<'PY'
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
coord={}
for l in open(f"{F}/am_dars.bed"):
    f=l.rstrip("\n").split("\t"); coord[f[3]]=(f[0],int(f[1]))
iv=collections.defaultdict(list); n=0
for i,l in enumerate(open(f"{F}/tfsc_hits.txt")):
    if not i or not l.strip(): continue
    p=l.rstrip("\n").split("\t")
    if len(p)<6 or p[2] not in coord: continue
    iv[(p[0],p[2])].append((int(p[3])-1,int(p[4]))); n+=1
o=open(f"{F}/tfsc_merged.bed","w"); k=0
for (m,e),v in iv.items():
    v.sort(); cur=list(v[0]); mg=[]
    for a,b in v[1:]:
        if a<=cur[1]: cur[1]=max(cur[1],b)
        else: mg.append(tuple(cur)); cur=[a,b]
    mg.append(tuple(cur)); c,s0=coord[e]
    for a,b in mg: k+=1; o.write(f"{c}\t{s0+a}\t{s0+b}\t{m}@{e}@{k}\n")
o.close(); print(f"  merged {k} intervals over {len(iv)} (motif,element) pairs")
PY
sort -k1,1 -k2,2n tfsc_merged.bed -o tfsc_merged.bed
for t in subtree222 wholetree; do
  bigWigAverageOverBed $H/${t}_perbase_lsd2.bw tfsc_merged.bed ${t}_tfsc_hits.tab
done
echo TFSCDONE
