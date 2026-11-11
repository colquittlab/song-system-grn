#!/usr/bin/env python
"""Can motif masking produce a neutral reference class inside enhancers?"""
import numpy as np, collections, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "mathtext.fontset":"custom","mathtext.rm":"Arial",
    "mathtext.it":"Arial:italic","mathtext.bf":"Arial:bold",
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2,
    "legend.handletextpad":.4,"legend.borderpad":.2,"legend.labelspacing":.3})
import matplotlib.pyplot as plt
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"; WARN="#7A6220"
TIERS=[("ar1",1,0.3),("sw37",37,8.4),("meta192",192,46.2),("all7446",7446,99.8)]
coord={l.rstrip("\n").split("\t")[3]:1 for l in open(f"{F}/loci_dars.bed")}
def load(f):
    d=collections.defaultdict(list)
    for l in open(f"{F}/{f}"):
        p=l.split("\t")
        if int(p[2])/max(int(p[1]),1)<.5: continue
        e=p[0].split("|")[0]
        if e in coord: d[e].append((float(p[5]),int(p[1])))
    return d
def wm(v): a=np.array(v); return float(np.average(a[:,0],weights=a[:,1]))
def ci(v,n=4000):
    b=np.array([np.mean(rng.choice(v,len(v))) for _ in range(n)]); return np.percentile(b,[2.5,97.5])
mot={}
for l in open(f"{F}/subtree222_motif2.tab"):
    p=l.split("\t")
    if int(p[2])/max(int(p[1]),1)>=.5: mot[p[0]]=float(p[5])
allc={l.split("\t")[3]:1 for l in open(f"{F}/allelem_cls.bed")}
h2e={}; n=0
for s in ("ra","g4"):
    for i,l in enumerate(open(f"{F}/fimo_motif2_{s}.txt")):
        if not i: continue
        e=l.split("\t")[2]
        if e not in allc: continue
        n+=1
        for loc in ("MAFB","EMX2","other"): h2e[f"motif2|{s}|{loc}|{n}"]=e
x=[t[2] for t in TIERS]
fig,axes=plt.subplots(1,2,figsize=(3.62,1.40))
fig.patch.set_facecolor("white")
ax=axes[0]
for tree,col,lab in (("wholetree","#8A929C","whole tree"),("subtree222","#1F3A5F","oscine subtree")):
    ys=[];los=[];his=[]
    for t,nm,pct in TIERS:
        D=load(f"{tree}_mcrest_{t}.tab"); v=np.array([wm(q) for q in D.values()])
        ys.append(v.mean()); lo,hi=ci(v); los.append(lo); his.append(hi)
    ax.fill_between(x,los,his,color=col,alpha=.15,linewidth=0)
    ax.plot(x,ys,color=col,lw=1.00,marker="o",ms=3.5,label=lab)
ax.axhline(0,color=INK,lw=0.50)
ax.axvspan(85,104,color=WARN,alpha=.10,linewidth=0,zorder=0)
ax.text(94.5,.44,"none\nleft",ha="center",va="top",fontsize=5,color=WARN,style="italic")
ax.set_xticks(x); ax.set_xticklabels([f"{p:.0f}%\n{n:,}" for _,n,p in TIERS],fontsize=5)
ax.set_xlabel("% masked  (motifs)",fontsize=6,color=INK)
ax.set_ylabel("phyloP, remainder",fontsize=6,color=INK)
ax.set_title("A   masking cannot reach neutrality",loc="left",fontsize=7,color=INK,
             fontweight="600",pad=10)
ax.set_ylim(-.02,.48); ax.set_xlim(-6,106)
ax.legend(frameon=False,fontsize=5,loc="lower left")
ax=axes[1]
for tag,col,lab in (("ra",RA,"Glut-RA"),("g4",G4,"GABA-4")):
    ys=[];los=[];his=[]
    for t,nm,pct in TIERS:
        S=load(f"subtree222_mcrest_{t}.tab")
        d=[v-wm(S[h2e[k]]) for k,v in mot.items()
           if k.split("|")[1]==tag and k.split("|")[2]=="MAFB"
           and h2e.get(k) in S]
        d=np.array(d); ys.append(d.mean()); lo,hi=ci(d); los.append(lo); his.append(hi)
    ax.fill_between(x,los,his,color=col,alpha=.15,linewidth=0)
    ax.plot(x,ys,color=col,lw=1.00,marker="o",ms=3.5,label=f"{lab} DARs at MAFB")
    dy=.055 if tag=="ra" else -.085
ax.axhline(0,color=INK,lw=0.50)
ax.axvspan(85,104,color=WARN,alpha=.10,linewidth=0,zorder=0)
ax.text(94.5,.62,"n.d.",ha="center",va="top",fontsize=5,color=WARN,style="italic")
ax.set_xlabel("% masked",fontsize=6,color=INK)
ax.set_ylabel("motif − remainder",fontsize=5,color=INK)
ax.set_title("B   the AR result strengthens as masking increases",loc="left",
             fontsize=7,color=INK,fontweight="600",pad=10)
ax.set_xlim(-6,106); ax.legend(frameon=False,fontsize=5,loc="upper left")
for a in axes:
    a.margins(x=.03,y=.05)
    a.locator_params(nbins=5)
    for sp in ("top","right"): a.spines[sp].set_visible(False)
    for sp in ("left","bottom"): a.spines[sp].set_color(RULE)
    a.tick_params(colors=INK3,labelsize=5)
fig.subplots_adjust(left=.085,right=.99,top=.885,bottom=.285,wspace=.16)
fig.savefig(f"{F}/results/compare_gw_mafb_emx2/fig_maskcurve.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb_emx2/fig_maskcurve.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_maskcurve.png")
