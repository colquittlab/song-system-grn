#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 16G
#SBATCH -c 2
#SBATCH -J sens
#SBATCH --array=0-21
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/sens/scan-%a.out
# One chunk of the motif panel through the same steps as s_tfsc.sh: FIMO over every
# Glut-RA/GABA-4 DAR, merge overlapping hits per (matrix, element), mean phyloP per
# merged interval on both trees.
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/ucsc_utilities:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus
H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload
C=$(printf "%02d" $SLURM_ARRAY_TASK_ID); cd $F
fimo --verbosity 1 --thresh 1e-4 --text sens/chunk$C.meme am_dars.fa > sens/chunk${C}_hits.txt
python3 - "$C" <<'PY'
import collections, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"; C=sys.argv[1]
coord={}
for l in open(f"{F}/am_dars.bed"):
    f=l.rstrip("\n").split("\t"); coord[f[3]]=(f[0],int(f[1]))
iv=collections.defaultdict(list)
for i,l in enumerate(open(f"{F}/sens/chunk{C}_hits.txt")):
    if not i or not l.strip(): continue
    p=l.rstrip("\n").split("\t")
    if len(p)<6 or p[2] not in coord: continue
    iv[(p[0],p[2])].append((int(p[3])-1,int(p[4])))
o=open(f"{F}/sens/chunk{C}_merged.bed","w"); k=0
for (m,e),v in iv.items():
    v.sort(); cur=list(v[0]); mg=[]
    for a,b in v[1:]:
        if a<=cur[1]: cur[1]=max(cur[1],b)
        else: mg.append(tuple(cur)); cur=[a,b]
    mg.append(tuple(cur)); c,s0=coord[e]
    for a,b in mg: k+=1; o.write(f"{c}\t{s0+a}\t{s0+b}\t{m}@{e}@{k}\n")
o.close(); print(f"  chunk {C}: {k} merged intervals over {len(iv)} (matrix,element) pairs")
PY
sort -k1,1 -k2,2n sens/chunk${C}_merged.bed -o sens/chunk${C}_merged.bed
for t in subtree222 wholetree; do
  bigWigAverageOverBed $H/${t}_perbase_lsd2.bw sens/chunk${C}_merged.bed sens/chunk${C}_${t}.tab
done
echo SENSDONE
