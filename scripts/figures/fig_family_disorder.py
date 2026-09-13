# Recovered from session records (the lab repository's reproducibility/recovered_scripts/family_disorder.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches. Per-family table -> results/tables/.
"""Group %IDR by 'Protein families' (first-token / superfamily level) across the 21,801 isoform set.
Rank families by median disorder; color by TF content."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root

df=pd.read_csv(REPO/"data/annotated/IDRisoforms_df_geneage.csv", low_memory=False)
df=df[df["Protein families"].notna()].copy()
df["family"]=df["Protein families"].astype(str).str.split(",").str[0].str.strip()
df["is_tf"]=df["tf_group"].astype(str)=="TF"

g=df.groupby("family").agg(
    n_iso=("pct_idr","size"),
    n_genes=("base_accession","nunique"),
    median_idr=("pct_idr","median"),
    mean_idr=("pct_idr","mean"),
    tf_frac=("is_tf","mean"),
).reset_index()

MIN_ISO=30
tab=g[g["n_iso"]>=MIN_ISO].sort_values("median_idr", ascending=False).reset_index(drop=True)
print(f"{len(tab)} families with >= {MIN_ISO} isoforms\n")
print("=== 15 MOST disordered families ===")
print(tab.head(15)[["family","n_iso","n_genes","median_idr","mean_idr","tf_frac"]].round(3).to_string(index=False))
print("\n=== 15 LEAST disordered families ===")
print(tab.tail(15)[["family","n_iso","n_genes","median_idr","mean_idr","tf_frac"]].round(3).to_string(index=False))

# save full table
outcsv=REPO/"results/tables/family_disorder_table.csv"; outcsv.parent.mkdir(parents=True, exist_ok=True)
g.sort_values("median_idr",ascending=False).round(3).to_csv(outcsv,index=False)
print("\nfull per-family table ->", outcsv)

# ---- figure: 16 most + 16 least disordered families (by median), extremes of the ranking ----
N_END=10
n_omit=max(0, len(tab)-2*N_END)
if len(tab)>2*N_END:
    show=pd.concat([tab.head(N_END), tab.tail(N_END)])
else:
    show=tab.copy()
show=show.reset_index(drop=True)

fig,ax=plt.subplots(figsize=(10.5, 0.42*len(show)+1.6))
ypos=np.arange(len(show))[::-1]   # highest median at top
cmap=matplotlib.colormaps["PuBu"]; norm=Normalize(0,1)
data=[df[df["family"]==f]["pct_idr"].dropna().values for f in show["family"]]
bp=ax.boxplot(data, positions=ypos, vert=False, widths=0.62,
              patch_artist=True, showfliers=False,
              medianprops=dict(color="black",linewidth=1.8))
for patch,tf in zip(bp["boxes"], show["tf_frac"]):
    patch.set_facecolor(cmap(norm(tf))); patch.set_edgecolor("black"); patch.set_linewidth(0.8)
import re
def short(f):
    s=re.sub(r'\s*\(TC[^)]*\)','',f)                       # drop transporter TC codes
    s=re.sub(r'\s+(super)?family$','',s,flags=re.I)        # drop trailing 'family'/'superfamily'
    s=s.replace("protein ","").replace("Krueppel C2H2-type zinc-finger","Krüppel C2H2 ZF")
    return s.strip()
ax.set_yticks(ypos)
ax.set_yticklabels([f"{short(f)}\n(n={n:,})" for f,n in zip(show["family"],show["n_iso"])],
                   fontsize=7.5, linespacing=0.9)
ax.tick_params(axis="y", pad=2, length=0)   # pull labels toward the plot
ax.set_xlabel("% IDR per isoform", fontsize=11)
ax.set_xlim(-2,102)
ax.axvline(df["pct_idr"].median(), color="gray", ls="--", lw=1, alpha=.7)
ax.text(df["pct_idr"].median()+1, len(show)-0.3, f"proteome median\n{df['pct_idr'].median():.0f}%",
        fontsize=7.5, color="gray", va="top")
# separator between the most- and least-disordered blocks
if n_omit>0:
    ybreak=ypos[N_END-1]-0.5
    ax.axhline(ybreak, color="gray", ls=":", lw=1)
ax.set_title(f"Intrinsic disorder varies markedly across protein families\n"
             f"({N_END} most vs {N_END} least disordered; families with ≥{MIN_ISO} isoforms)",
             fontsize=12.5, fontweight="bold", pad=10)
sm=cm.ScalarMappable(cmap=cmap,norm=norm); sm.set_array([])
cb=fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.03)
cb.set_label("fraction of isoforms that are TFs", fontsize=9)
fig.subplots_adjust(left=0.16, right=0.99, top=0.94, bottom=0.08)
fig.savefig(REPO/"figures/proteome/fig_family_disorder.pdf", bbox_inches="tight")
print("figure saved.")
