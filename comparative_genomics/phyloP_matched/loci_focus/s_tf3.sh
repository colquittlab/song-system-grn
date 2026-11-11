#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 48G
#SBATCH -c 8
#SBATCH -J tf3
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/tf3-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/ucsc_utilities:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; cd $F
cat tf_Fosl2.meme > tf3.meme
for m in Mafb Nfia; do sed -n '/^MOTIF/,$p' tf_$m.meme >> tf3.meme; done
echo "  motifs: $(grep -c '^MOTIF' tf3.meme)"
[ -s tf3_hits.txt ] || fimo --verbosity 1 --thresh 1e-4 --text tf3.meme am_dars.fa > tf3_hits.txt
echo "  hits: $(( $(wc -l < tf3_hits.txt) - 1 ))"
python3 - <<'PY'
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
coord={}
for l in open(f"{F}/am_dars.bed"):
    f=l.rstrip("\n").split("\t"); coord[f[3]]=(f[0],int(f[1]))
iv=collections.defaultdict(list); n=0
for i,l in enumerate(open(f"{F}/tf3_hits.txt")):
    if not i or not l.strip(): continue
    p=l.rstrip("\n").split("\t")
    if len(p)<6: continue
    if p[2] not in coord: continue
    iv[(p[0],p[2])].append((int(p[3])-1,int(p[4]))); n+=1
o=open(f"{F}/tf3_merged.bed","w"); k=0
for (m,e),v in iv.items():
    v.sort(); cur=list(v[0]); mg=[]
    for a,b in v[1:]:
        if a<=cur[1]: cur[1]=max(cur[1],b)
        else: mg.append(tuple(cur)); cur=[a,b]
    mg.append(tuple(cur)); c,s0=coord[e]
    for a,b in mg: k+=1; o.write(f"{c}\t{s0+a}\t{s0+b}\t{m}@{e}@{k}\n")
o.close(); print(f"  merged: {k} over {len(iv)} pairs")
PY
sort -k1,1 -k2,2n tf3_merged.bed -o tf3_merged.bed
for t in subtree222 wholetree; do
  bigWigAverageOverBed $H/${t}_perbase_lsd2.bw tf3_merged.bed ${t}_tf3_hits.tab
done
echo TF3DONE
