# Recovered from session records (the lab repository's reproducibility/recovered_scripts/fig15_tf_nontf_remake_inline.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches. Output renamed (see below).
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]  # repository root
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

d=pd.read_csv(REPO/'data/genomic/analysis/idr_gain_loss_by_tf.csv').set_index('cohort')
rows=['TF','NON-TF']
CONS="#2CA089"; REM="#2E7EBB"; ADD="#E0A32E"
cats=[('conserved',CONS,'conserved (in both)'),('removed',REM,'removed (canonical only)'),('added',ADD,'added (alternative only)')]

# chi-square TF vs non-TF fate
tab=d.loc[rows,['conserved','removed','added']].values
chi2,pval,_,_=chi2_contingency(tab)
n=tab.sum(); cram=np.sqrt(chi2/(n*(min(tab.shape)-1)))
print(f"TF vs non-TF fate: chi2={chi2:.2f} P={pval:.2f} Cramer V={cram:.3f}")

fig,(axA,axB)=plt.subplots(1,2,figsize=(12,4),gridspec_kw={'width_ratios':[1.25,1]})
# ---- Panel A: stacked % horizontal ----
y=np.arange(len(rows))[::-1]
left=np.zeros(len(rows))
for key,col,lab in cats:
    vals=d.loc[rows,f'pct_{key}'].values
    axA.barh(y,vals,left=left,color=col,edgecolor='white',height=0.62,label=lab)
    for yi,(v,l) in enumerate(zip(vals,left)):
        axA.text(l+v/2, y[yi], f"{v:.0f}%", ha='center',va='center',
                 color='white',fontweight='bold',fontsize=11)
    left+=vals
axA.set_yticks(y); axA.set_yticklabels(rows,fontsize=12,fontweight='bold')
axA.set_xlim(0,100); axA.set_xlabel("% of IDR exons (per canonical↔alternative comparison)",fontsize=10)
axA.set_title("a",loc='left',fontsize=14,fontweight='bold')
axA.legend(ncol=3,frameon=False,fontsize=9.5,loc='upper center',bbox_to_anchor=(0.5,1.16))
axA.spines[['top','right']].set_visible(False)

# ---- Panel B: grouped counts, log ----
x=np.arange(len(rows)); w=0.26
for i,(key,col,lab) in enumerate(cats):
    vals=d.loc[rows,key].values
    bars=axB.bar(x+(i-1)*w,vals,w,color=col,label=lab.split(' (')[0])
    for b,v in zip(bars,vals):
        axB.text(b.get_x()+b.get_width()/2,v*1.08,f"{v:,}",ha='center',fontsize=8.5)
axB.set_yscale('log'); axB.set_ylim(100,1e5)
axB.set_xticks(x); axB.set_xticklabels(rows,fontsize=12,fontweight='bold')
axB.set_ylabel("number of IDR exons (log scale)",fontsize=10)
axB.set_title("b",loc='left',fontsize=14,fontweight='bold')
axB.legend(frameon=False,fontsize=9,loc='upper right')
axB.spines[['top','right']].set_visible(False)

fig.suptitle("Alternative splicing leaves most IDR exons intact; removal exceeds addition — equally in TFs and non-TFs",
             fontsize=13,fontweight='bold',y=1.04)
fig.text(0.5,-0.04,f"TF vs non-TF fate distribution: χ² test on 2×3 (TF/non-TF × conserved/removed/added): "
         f"P = {pval:.2f} (n.s.), Cramér's V = {cram:.3f} (negligible).",
         ha='center',fontsize=8.5,color='0.35')
fig.tight_layout()
# RENAMED OUTPUT: this TF/non-TF-only version used to overwrite s10b's fig15 PDF (while s10b's PNG
# kept the ALL/TF/non-TF version). It now has its own name so s10b's outputs stay consistent.
fig.savefig(REPO/'figures/genomic/fig15_idr_gain_loss_tf_vs_nontf.pdf',bbox_inches='tight')
print("saved fig15_idr_gain_loss_tf_vs_nontf.pdf")