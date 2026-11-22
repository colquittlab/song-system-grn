#!/usr/bin/env python
"""Per-instance phyloP across the JASPAR AR motif: every site as its own row.

The position-averaged profile can only show the mean shape; this asks whether
that mean is a property of most sites or of a constrained subset.
"""
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
from matplotlib.gridspec import GridSpec
import pandas as pd, logomaker
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
W=16; VM=2.0
txt=open(f"{F}/ar_jaspar.meme").read()
blk=txt.split("letter-probability")[1].split("\n")[1:]
pwm=np.array([[float(x) for x in l.split()] for l in blk if len(l.split())==4][:W])
ic=np.array([2+(p*np.log2(np.clip(p,1e-9,None))).sum() for p in pwm])
info=pd.DataFrame(pwm*ic[:,None],columns=list("ACGT"))
inst=collections.defaultdict(dict)
for l in open(f"{F}/subtree222_jaspos.tab"):
    p=l.split("\t")
    if int(p[2])==0: continue
    s,loc,k,idx=p[0].split("|"); inst[(s,loc,idx)][int(k)]=float(p[5])
def M(s,loc=None):
    V=[[v[k] for k in range(W)] for (ss,ll,ix),v in inst.items()
       if ss==s and (loc is None or ll==loc) and len(v)==W]
    V=np.array(V)
    return V[np.argsort(-V.mean(1))]          # most-constrained sites at the top
PAN=[("Glut-RA DARs, genome-wide","ra",None,RA),
     ("GABA-4 DARs, genome-wide","g4",None,G4),
     ("Glut-RA DARs at MAFB","ra","MAFB",RA),
     ("GABA-4 DARs at MAFB","g4","MAFB",G4)]
mats=[(t,M(s,l),c) for t,s,l,c in PAN]
hr=[1.05]+[max(0.9,np.log10(len(m)+1)*2.0) for _,m,_ in mats]
fig=plt.figure(figsize=(2.48,3.20)); fig.patch.set_facecolor("white")
gs=GridSpec(len(hr),2,height_ratios=hr,width_ratios=[26,1],hspace=.80,wspace=.16,
            left=.145,right=.90,top=.958,bottom=.135)
ax=fig.add_subplot(gs[0,0])
logomaker.Logo(info,ax=ax,color_scheme="classic",show_spines=False,font_name="Arial")
ax.set_ylabel("bits",fontsize=6,color=INK); ax.set_yticks([0,2])
ax.set_xticks([]); ax.tick_params(colors=INK3,labelsize=5)
ax.set_title("JASPAR MA0007.3 (AR)",loc="left",fontsize=5,color=INK,fontweight="600",pad=6)
im=None
for i,(t,m,c) in enumerate(mats):
    ax=fig.add_subplot(gs[i+1,0])
    im=ax.imshow(m,aspect="auto",cmap="RdBu_r",vmin=-VM,vmax=VM,interpolation="nearest")
    ax.set_ylabel(f"{len(m)} sites",fontsize=5,color=INK3)
    ax.set_yticks([])
    ax.set_xticks(range(W))
    ax.set_xticklabels(range(1,W+1) if i==len(mats)-1 else [],fontsize=5,color=INK3)
    ax.set_title(t,loc="left",fontsize=7,color=c,fontweight="600",pad=5)
    for sp in ax.spines.values(): sp.set_color(RULE)
    # marginal: per-site mean, to show the spread the heatmap rows are sorted on
    axm=fig.add_subplot(gs[i+1,1])
    axm.plot(m.mean(1),np.arange(len(m)),color=c,lw=0.55)
    axm.axvline(0,color=RULE,lw=.8)
    axm.set_ylim(len(m)-.5,-.5); axm.set_xlim(-VM,VM)
    axm.set_xticks([]); axm.set_yticks([])
    for sp in axm.spines.values(): sp.set_color(RULE)
    if i==0: axm.set_title("site\nmean",fontsize=7,color=INK3,pad=5)
ax.set_xlabel("motif position",fontsize=6,color=INK,labelpad=6)
cax=fig.add_axes([.145,.062,.28,.010])
cb=fig.colorbar(im,cax=cax,orientation="horizontal")
cb.set_label("oscine phyloP",fontsize=5,color=INK,labelpad=3)
cb.ax.tick_params(colors=INK3,labelsize=5)
cb.outline.set_edgecolor(RULE)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_heatmap.png",dpi=190,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_heatmap.pdf",dpi=190,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar_heatmap.png")
for t,m,_ in mats:
    r=m.mean(1)
    print(f"  {t:<32} n={len(m):<5} site means: median {np.median(r):+.3f}  "
          f"IQR [{np.percentile(r,25):+.3f},{np.percentile(r,75):+.3f}]  "
          f"frac>0.5 {np.mean(r>0.5):.2f}  frac<0 {np.mean(r<0):.2f}")
