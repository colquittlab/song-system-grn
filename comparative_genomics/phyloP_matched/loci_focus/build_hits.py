"""Merge FIMO hits within each (motif, element) so overlapping calls count once."""
import collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
coord={}
for l in open(f"{F}/am_dars.bed"):
    f=l.rstrip("\n").split("\t"); coord[f[3]]=(f[0],int(f[1]))
iv=collections.defaultdict(list)
n=0
for i,l in enumerate(open(f"{F}/am_hits.txt")):
    if not i or not l.strip(): continue
    p=l.rstrip("\n").split("\t")
    if len(p)<6: continue
    m,e=p[0],p[2]
    if e not in coord: continue
    iv[(m,e)].append((int(p[3])-1,int(p[4]))); n+=1
print(f"  parsed {n} hits over {len(iv)} (motif,element) pairs",flush=True)
out=open(f"{F}/am_merged.bed","w"); k=0
for (m,e),v in iv.items():
    v.sort(); cur=list(v[0]); mg=[]
    for a,b in v[1:]:
        if a<=cur[1]: cur[1]=max(cur[1],b)
        else: mg.append(tuple(cur)); cur=[a,b]
    mg.append(tuple(cur))
    c,s0=coord[e]
    for a,b in mg:
        k+=1; out.write(f"{c}\t{s0+a}\t{s0+b}\t{m}@{e}@{k}\n")
out.close()
print(f"  merged intervals: {k}",flush=True)
