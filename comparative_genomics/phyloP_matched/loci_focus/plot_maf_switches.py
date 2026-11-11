#!/usr/bin/env python
"""Per-site oscine/outgroup base switches across AR motif sites at MAFB.

One row per site, one cell per motif column. Cell shade is how much the oscine
major allele gains in frequency relative to the outgroup; cells where the two
groups take a different majority base, with >=80% support on each side, are
labelled with the change.
"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2})
import matplotlib.pyplot as plt
import pandas as pd, logomaker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
def homer_pwm(f): return np.array([[float(x) for x in l.split()] for l in open(f)
                                   if not l.startswith(">")])
def meme_pwm(f,W):
    b=open(f).read().split("letter-probability")[1].split("\n")[1:]
    return np.array([[float(x) for x in l.split()] for l in b if len(l.split())==4][:W])
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
MINO,MING,THR=6,10,0.80
CMAP=LinearSegmentedColormap.from_list("sh",["#FFFFFF","#FBE3C8","#F2A25C","#C7521B"])
R=json.load(open(f"{F}/results/mafb_1p5mb/maf_ar_mafb.json"))
cls=lambda s:"Glut-RA" if s.startswith("RA_") else "GABA-4"
short=lambda s:s.split("_",1)[1].rsplit(".",1)[-1].replace("_",".")
def cellvals(r):
    d=[];lab=[];ok=[]
    for c in r["cols"]:
        good=c["oscN"]>=MINO and c["outN"]>=MING
        ok.append(good)
        if not good: d.append(np.nan); lab.append(""); continue
        # frequency of the oscine major allele, in each group
        sw = c["osc"]!=c["out"] and c["oscF"]>=THR and c["outF"]>=THR
        a=c["osc"]                                  # the oscine major allele
        fo=c["oscC"].get(a,0)/max(c["oscN"],1)
        fg=c["outC"].get(a,0)/max(c["outN"],1)
        d.append(max(0.0,fo-fg))                    # one definition for every cell
        lab.append(f"{c['out']}→{c['osc']}" if sw else "")
    return np.array(d,float),lab,ok
def nsw(r): return sum(1 for c in r["cols"]
                       if c["oscN"]>=MINO and c["outN"]>=MING and c["osc"]!=c["out"]
                       and c["oscF"]>=THR and c["outF"]>=THR)
PAN=[("full","ARE-full",12),("jaspar","JASPAR AR",16)]
PWM={"full":homer_pwm(f"{HM}/motif2.motif"),"jaspar":meme_pwm(f"{F}/ar_jaspar.meme",16)}
fig=plt.figure(figsize=(6.8,4.45)); fig.patch.set_facecolor("white")
gs=fig.add_gridspec(2,2,height_ratios=[1,13],width_ratios=[12,16],hspace=.045,wspace=.30)
for ci,(mk,ml,W) in enumerate(PAN):
    axl=fig.add_subplot(gs[0,ci])
    ic=np.array([2+(q*np.log2(np.clip(q,1e-9,None))).sum() for q in PWM[mk]])
    logomaker.Logo(pd.DataFrame(PWM[mk]*ic[:,None],columns=list("ACGT")),
                   ax=axl,color_scheme="classic",show_spines=False,font_name="Arial")
    axl.set_xlim(-.5,W-.5); axl.set_ylim(0,2); axl.set_xticks([]); axl.set_yticks([0,2])
    axl.tick_params(colors=INK3,labelsize=4); axl.set_ylabel("bits",fontsize=5,color=INK)
    ax=fig.add_subplot(gs[1,ci])
    S=[r for r in R if r["motif"]==mk]
    S.sort(key=lambda r:(cls(r["site"])!="Glut-RA",-nsw(r),short(r["site"])))
    M=np.full((len(S),W),np.nan); L=[]
    for i,r in enumerate(S):
        d,lab,ok=cellvals(r); M[i,:len(d)]=d; L.append(lab)
    ax.imshow(M,cmap=CMAP,vmin=0,vmax=1,aspect="auto",interpolation="nearest",
              extent=(-.5,W-.5,len(S)-.5,-.5))
    ax.set_xlim(-.5,W-.5)
    for i in range(len(S)):
        for j in range(W):
            if np.isnan(M[i,j]):
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor="#EDEFF2",
                                       edgecolor="none",zorder=2))
            elif L[i][j]:
                ax.text(j,i,L[i][j],ha="center",va="center",fontsize=3.4,
                        color="#3A1A06",fontweight="700",zorder=3,clip_on=False)
    ax.set_xticks(range(W)); ax.set_xticklabels(range(1,W+1),fontsize=4.4)
    ax.set_yticks(range(len(S)))
    ax.set_yticklabels([short(r["site"]) for r in S],fontsize=3.9)
    ax.set_xlabel("motif position",fontsize=6,color=INK)
    axl.set_title(f"{'AB'[ci]}   {ml} ({W} bp), {len(S)} sites at MAFB",
                  loc="left",fontsize=7,color=INK,fontweight="600",pad=2)
    nRA=sum(1 for r in S if cls(r["site"])=="Glut-RA")
    ax.axhline(nRA-.5,color=INK,lw=.7)
    ax.set_xticks(np.arange(-.5,W,1),minor=True)
    ax.set_yticks(np.arange(-.5,len(S),1),minor=True)
    ax.grid(which="minor",color="white",lw=.45); ax.tick_params(which="minor",length=0)
    ax.tick_params(colors=INK3,labelsize=4.4)
    for lb,r in zip(ax.get_yticklabels(),S):            # after tick_params, which resets colour
        lb.set_color(RA if cls(r["site"])=="Glut-RA" else G4); lb.set_fontsize(3.9)
    for sp in ax.spines.values(): sp.set_visible(False)
sm=plt.cm.ScalarMappable(cmap=CMAP,norm=plt.Normalize(0,1))
cb=fig.colorbar(sm,ax=fig.axes,fraction=.014,pad=.028,aspect=30)
cb.set_label("frequency of the oscine major allele:\n13 oscines minus 26 outgroup",fontsize=5,color=INK)
cb.ax.tick_params(labelsize=4.4,colors=INK3,width=.4,length=1.5)
cb.outline.set_visible(False)
fig.legend(handles=[Line2D([],[],marker="s",color=RA,ls="",ms=3,label="Glut-RA DAR"),
                    Line2D([],[],marker="s",color=G4,ls="",ms=3,label="GABA-4 DAR"),
                    Line2D([],[],marker="s",color="#EDEFF2",ls="",ms=3,
                           label="too few aligned genomes")],
           frameon=False,fontsize=4.6,loc="lower left",ncol=3,bbox_to_anchor=(.07,-.008))
fig.subplots_adjust(left=.10,right=.885,top=.945,bottom=.075)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/mafb_1p5mb/fig_maf_switches.{e}",dpi=220,facecolor="white",
                bbox_inches="tight",pad_inches=.015)
print("  wrote fig_maf_switches")
