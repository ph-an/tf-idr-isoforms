# Recovered from session records (the lab repository's reproducibility/recovered_scripts/fig16b_idrfrac_by_age.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches.
"""fig16b -- fig16's CONTINUOUS per-exon idr_frac, recast as split violins across the four
gene-age bins, Non-TF (blue) vs TF (orange). Alternative exons, unique per gene+exon_id
(same dedup as fig16). Solid=median, dashed=quartiles (fig16 style); overlaid mean markers +
trend line make the evolutionary shift readable where the bimodal medians saturate at 0/1."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, gaussian_kde
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root

# ---- fig16 source: unique ALTERNATIVE exons, continuous idr_frac ----
ex = pd.read_csv(REPO/"data/genomic/analysis/exon_idr_annotation.csv")
ex = ex[ex.usage.isin(["alternative","constitutive"])]
u = (ex.groupby(["gene","exon_id"])
       .agg(idr_frac=("idr_frac","mean"), usage=("usage","first"), is_tf=("is_tf","first")).reset_index())
alt = u[u.usage=="alternative"].copy()
m = pd.read_csv(REPO/"data/annotated/IDRisoforms_df_geneage.csv",
    usecols=["base_accession","age_bin","tf_group"], low_memory=False).drop_duplicates("base_accession")
df = alt.merge(m, left_on="gene", right_on="base_accession", how="left").dropna(subset=["age_bin","tf_group"])

def stars(p):
    if pd.isna(p): return "NA"
    return "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "ns"
def bracket(ax,x1,x2,y,text,color="black",tick=0.018,lw=1.2,fs=12):
    ax.plot([x1,x1,x2,x2],[y-tick,y,y,y-tick],color=color,lw=lw,clip_on=False,solid_capstyle="butt")
    ax.text((x1+x2)/2,y+0.006,text,ha="center",va="bottom",color=color,fontsize=fs,fontweight="bold",clip_on=False)
def darken(hexc,f=0.55):
    r,g,b=mc.to_rgb(hexc); return (r*f,g*f,b*f)

age_order=["<100","100-500","500-1000",">1000"]
df["age_bin"]=pd.Categorical(df["age_bin"],categories=age_order,ordered=True)

INK,MUTED,GRID="#1a1a1a","#666","#E6E6E6"
plt.rcParams.update({"font.family":"sans-serif","font.size":10,"axes.titlesize":11.5,
    "axes.titleweight":"bold","axes.edgecolor":"#444","axes.linewidth":0.9,"text.color":INK,
    "xtick.color":"#444","ytick.color":"#444","axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(figsize=(9.6,5.6))
ax.grid(axis="y",color=GRID,lw=.9); ax.set_axisbelow(True)
centers=np.arange(1,len(age_order)+1); xgap=0.06
non_pos=centers-xgap/2; tf_pos=centers+xgap/2
vw=0.82; mhw=vw/2
NONTF="#4C72B0"; TF="#DD8452"; NONTFd=darken(NONTF); TFd=darken(TF)

def half_violin(arr,pos,side,face,dark):
    arr=np.asarray(arr,float)
    if len(arr)<2 or np.std(arr)==0: return
    vp=ax.violinplot([arr],positions=[pos],widths=vw,showextrema=False,
                     side="low" if side=="left" else "high")
    for b in vp["bodies"]:
        b.set_facecolor(face); b.set_edgecolor("black"); b.set_linewidth(1.0); b.set_alpha(0.72)
        v=b.get_paths()[0].vertices; v[:,1]=np.clip(v[:,1],0,1)
    q1,med,q3=np.percentile(arr,[25,50,75])
    kde=gaussian_kde(arr); maxd=kde(np.linspace(arr.min(),arr.max(),300)).max()
    wq1,wmed,wq3=(mhw*(d/maxd) for d in kde([q1,med,q3]))
    for yv,wd in [(q1,wq1),(q3,wq3)]:
        lo,hi=(pos-wd,pos) if side=="left" else (pos,pos+wd)
        ax.hlines(yv,lo,hi,lw=1.5,color=dark,linestyles=(0,(4,2)),zorder=5)
    lo,hi=(pos-wmed,pos) if side=="left" else (pos,pos+wmed)
    ax.hlines(med,lo,hi,lw=3.0,color=dark,zorder=6)

pvals={}; ticklabels=[]; non_means=[]; tf_means=[]
for i,age in enumerate(age_order):
    nv=df.query("age_bin==@age and tf_group=='Non-TF'")["idr_frac"].to_numpy(float)
    tv=df.query("age_bin==@age and tf_group=='TF'")["idr_frac"].to_numpy(float)
    half_violin(nv,non_pos[i],"left",NONTF,NONTFd)
    half_violin(tv,tf_pos[i],"right",TF,TFd)
    non_means.append(nv.mean()); tf_means.append(tv.mean())
    pvals[age]=mannwhitneyu(nv,tv,alternative="two-sided")[1] if len(nv) and len(tv) else np.nan
    ticklabels.append(f"{age} Ma\nNon-TF={len(nv):,}\nTF={len(tv):,}")

# mean trend lines + markers (readable where bimodal medians saturate)
ax.plot(non_pos,non_means,"-",color=NONTFd,lw=1.8,zorder=7,alpha=0.9)
ax.plot(tf_pos,tf_means,"-",color=TFd,lw=1.8,zorder=7,alpha=0.9)
ax.scatter(non_pos,non_means,s=52,facecolor="white",edgecolor=NONTFd,lw=1.8,zorder=8)
ax.scatter(tf_pos,tf_means,s=52,facecolor="white",edgecolor=TFd,lw=1.8,zorder=8)

ax.axhline(0.5,ls=(0,(6,4)),lw=1.2,color="#777",zorder=1)

ax.set_xticks(centers); ax.set_xticklabels(ticklabels,fontsize=8.5)
ax.set_xlabel("Gene age bin",labelpad=8)
ax.set_ylabel("per-exon IDR fraction\n(continuous idr_frac, alternative exons)")
ax.set_xlim(0.4,len(age_order)+0.6); ax.set_ylim(-0.05,1.16)
ax.set_yticks([0,0.25,0.5,0.75,1.0])
ax.set_title("Alternative-exon disorder diverges over evolution: TFs rise, non-TFs fall (fig16 metric)")
for i,age in enumerate(age_order):
    bracket(ax,centers[i]-0.24,centers[i]+0.24,1.075,stars(pvals[age]))

lg=[Patch(facecolor=NONTF,edgecolor="black",alpha=0.72,label="Non-TF"),
    Patch(facecolor=TF,edgecolor="black",alpha=0.72,label="TF"),
    Line2D([0],[0],color="#555",lw=1.8,marker="o",markerfacecolor="white",
           markeredgecolor="#555",markersize=7,label="group mean (trend)")]
ax.legend(handles=lg,frameon=True,edgecolor="black",facecolor="white",framealpha=1.0,
          loc="upper left",bbox_to_anchor=(1.004,1.0),borderaxespad=0,fontsize=8.5,handlelength=1.4)
fig.text(0.5,0.006,"Continuous per-exon idr_frac (fig16), unique alternative exons (dedup gene+exon_id). "
         "Solid=median, dashed=quartiles, grey dashed line = IDR-exon cutoff (0.5); distribution is bimodal (0 or 1) "
         "so means/trend lines carry the shift.  ns p\u22650.05   ****p<0.0001",
         ha="center",va="bottom",fontsize=7.1,color=MUTED)
fig.subplots_adjust(left=0.115,right=0.985,top=0.90,bottom=0.24)

OUT=REPO/"figures/genomic"
fig.savefig(OUT/"fig16b_idrfrac_by_age.pdf",bbox_inches="tight")
print("Non-TF means:", [round(x,2) for x in non_means])
print("TF means:    ", [round(x,2) for x in tf_means])
print("saved fig16b_idrfrac_by_age.pdf")
