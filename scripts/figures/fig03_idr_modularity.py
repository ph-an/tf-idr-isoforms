# Recovered from session records (the lab repository's reproducibility/recovered_scripts/fig03_idr_modularity_inline.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches.
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]  # repository root
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

idr=pd.read_parquet(REPO/'data/genomic/interim/idr_exon_overlap_table.parquet')
non=pd.read_parquet(REPO/'data/genomic/interim/nonidr_exon_overlap_table.parquet')
col='n_overlapping_coding_exons'
def dist(df):
    n=df[col].dropna()
    return np.array([(n==1).mean(), (n==2).mean(), (n>=3).mean()])*100
di=dist(idr); dn=dist(non)
print("IDR  1/2/>=3 %:", di.round(1), "n=",len(idr))
print("nonIDR 1/2/>=3 %:", dn.round(1), "n=",len(non))

IDR_ORANGE="#E4761B"; NON_BLUE="#4C72B0"
labels=["1 exon","2 exons","≥3 exons"]; xpos=np.arange(3); w=0.38
fig,ax=plt.subplots(figsize=(8,5.2))
b1=ax.bar(xpos-w/2, di, w, color=IDR_ORANGE, label="IDR segments", edgecolor="none")
b2=ax.bar(xpos+w/2, dn, w, color=NON_BLUE, label="non-IDR (ordered)", edgecolor="none")
for bars,vals in [(b1,di),(b2,dn)]:
    for r,v in zip(bars,vals):
        ax.text(r.get_x()+r.get_width()/2, v+1.2, f"{v:.0f}%", ha="center", fontsize=11)
ax.set_xticks(xpos); ax.set_xticklabels(labels, fontsize=12)
ax.set_ylabel("% of segments", fontsize=12)
ax.set_ylim(0,90)
ax.set_title("IDR modularity: exons spanned per segment", fontsize=15, fontweight="bold", loc="left")
ax.legend(frameon=False, fontsize=12, loc="upper left", bbox_to_anchor=(0.02,0.98))
ax.spines[['top','right']].set_visible(False)
ax.grid(axis='y', color="0.9", lw=0.8); ax.set_axisbelow(True)
fold=di[0]/dn[0]
fig.text(0.5,-0.01, f"IDRs are ~{fold:.0f}× more likely than ordered regions to be encoded by a single exon",
         ha="center", fontsize=10.5, style="italic", color="0.3")
fig.tight_layout()
for ext in ["pdf"]:
    fig.savefig(REPO/f"figures/genomic/fig03_idr_modularity.{ext}", bbox_inches="tight")
print("saved fig03_idr_modularity.pdf")