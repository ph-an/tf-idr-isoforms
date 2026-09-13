# Recovered from session records (the lab repository's reproducibility/recovered_scripts/fig09_violin.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches;
# alternative exons = exon_usage == "alternative" (mapped transcripts, the paper definition) instead of
# GENCODE-wide is_constitutive == False (2026-09-13; old version in the lab repository's archive/figures_superseded/).
"""fig09 as a fig4a-style split violin: per-gene % of alternative exons encoding an IDR,
by gene age, Non-TF (blue) vs TF (orange). Reveals opposing trends the pooled bars hid."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, gaussian_kde
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root

# ---- per-gene alt-exon IDR fraction ----
ex=pd.read_csv(REPO/"data/genomic/processed/exon_level_table.csv",
   usecols=["base_accession","exon_usage","is_idr_exon","is_coding_exon"], low_memory=False)
alt=ex[(ex["is_coding_exon"]==True)&(ex["exon_usage"]=="alternative")]
g=alt.groupby("base_accession").agg(frac=("is_idr_exon","mean"), n_alt=("is_idr_exon","size")).reset_index()
g["frac"]*=100
g=g[g["n_alt"]>=3]
m=pd.read_csv(REPO/"data/annotated/IDRisoforms_df_geneage.csv",
   usecols=["base_accession","age_bin","tf_group"], low_memory=False).drop_duplicates("base_accession")
df=g.merge(m,on="base_accession",how="left").dropna(subset=["age_bin","tf_group"])

def format_p(p):
    return "p < 1e-300" if p<1e-300 else (f"p = {p:.2e}" if p<0.001 else f"p = {p:.3f}")
def stars(p):
    if pd.isna(p): return "NA"
    return "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "ns"
def bracket(ax,x1,x2,y,text,color="black",tick=2.0,lw=1.2,fs=12):
    ax.plot([x1,x1,x2,x2],[y-tick,y,y,y-tick],color=color,lw=lw,clip_on=False,solid_capstyle="butt")
    ax.text((x1+x2)/2,y+0.6,text,ha="center",va="bottom",color=color,fontsize=fs,fontweight="bold",clip_on=False)

age_order=["<100","100-500","500-1000",">1000"]
metric="frac"
plot_df=df.copy()
plot_df["age_bin"]=pd.Categorical(plot_df["age_bin"],categories=age_order,ordered=True)

fig,ax=plt.subplots(figsize=(9,5.2))     # matches fig4a dimensions
centers=np.arange(1,len(age_order)+1); xgap=0.05
non_pos=centers-xgap/2; tf_pos=centers+xgap/2
vw=0.75; mhw=vw/2
NONTF="#4C72B0"; TF="#DD8452"; NONTFd="#2F4B7C"; TFd="#A95A2C"
pvals={}
def kde_w(vals,y):
    vals=np.asarray(pd.Series(vals).dropna())
    if len(vals)<2 or np.std(vals)==0: return mhw*0.2
    k=gaussian_kde(vals); yg=np.linspace(vals.min(),vals.max(),300)
    return mhw*(k([y])[0]/k(yg).max())
def qmed(x,vals,side,dark):
    vals=pd.Series(vals).dropna()
    if len(vals)==0: return
    q1,med,q3=np.percentile(vals,[25,50,75])
    for yv in [q1,q3]:
        w=kde_w(vals,yv); ax.hlines(yv,x-w if side=="left" else x,x if side=="left" else x+w,lw=1.1,color=dark,ls="dashed")
    w=kde_w(vals,med); ax.hlines(med,x-w if side=="left" else x,x if side=="left" else x+w,lw=2.4,color=dark)

ticklabels=[]
for i,age in enumerate(age_order):
    nv=plot_df.query("age_bin==@age and tf_group=='Non-TF'")[metric].dropna().to_numpy(dtype=float)
    tv=plot_df.query("age_bin==@age and tf_group=='TF'")[metric].dropna().to_numpy(dtype=float)
    if len(nv)>=2 and np.std(nv)>0:
        for b in ax.violinplot([nv],positions=[non_pos[i]],widths=vw,showextrema=False,side="low")["bodies"]:
            b.set_facecolor(NONTF); b.set_edgecolor("black"); b.set_linewidth(1.0); b.set_alpha(0.75)
        qmed(non_pos[i],nv,"left",NONTFd)
    if len(tv)>=2 and np.std(tv)>0:
        for b in ax.violinplot([tv],positions=[tf_pos[i]],widths=vw,showextrema=False,side="high")["bodies"]:
            b.set_facecolor(TF); b.set_edgecolor("black"); b.set_linewidth(1.0); b.set_alpha(0.75)
        qmed(tf_pos[i],tv,"right",TFd)
    pvals[age]=mannwhitneyu(nv,tv,alternative="two-sided")[1] if len(nv) and len(tv) else np.nan
    ticklabels.append(f"{age} Ma\nNon-TF={len(nv):,}\nTF={len(tv):,}")

ax.set_xticks(centers); ax.set_xticklabels(ticklabels,fontsize=8.5)
ax.set_xlabel("Gene age bin",labelpad=8); ax.set_ylabel("% of alt. exons encoding an IDR (per gene)")
ax.set_xlim(0.4,len(age_order)+0.6); ax.set_ylim(-4,122)
ax.set_title("IDR-targeted splicing by gene age: opposing trends in TFs vs Non-TFs")
for i,age in enumerate(age_order):
    bracket(ax,centers[i]-0.24,centers[i]+0.24,104,stars(pvals[age]))
lg=[Patch(facecolor=NONTF,edgecolor="black",alpha=0.75,label="Non-TF"),
    Patch(facecolor=TF,edgecolor="black",alpha=0.75,label="TF")]
ax.legend(handles=lg,frameon=True,edgecolor="black",facecolor="white",framealpha=1.0,
          loc="upper right",bbox_to_anchor=(0.988,0.985),borderaxespad=0,fontsize=9,ncol=2,columnspacing=1.0,handlelength=1.3)
fig.text(0.5,0.005,"ns  p≥0.05      *  p<0.05      **  p<0.01      ***  p<0.001      ****  p<0.0001",
         ha="center",va="bottom",fontsize=7.5,color="#444444")
fig.subplots_adjust(left=0.10,right=0.985,top=0.91,bottom=0.24)
for ext in ["pdf"]:
    fig.savefig(REPO/f"figures/genomic/fig09_gene_age.{ext}",bbox_inches="tight")
print("medians Non-TF:", plot_df[plot_df.tf_group=='Non-TF'].groupby('age_bin')[metric].median().round(1).to_dict())
print("medians TF:    ", plot_df[plot_df.tf_group=='TF'].groupby('age_bin')[metric].median().round(1).to_dict())
print("saved (overwrote fig09 pdf)")
