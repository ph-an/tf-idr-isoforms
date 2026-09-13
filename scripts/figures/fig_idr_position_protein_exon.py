# Recovered from session records (the lab repository's reproducibility/recovered_scripts/idr_position.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches.
# DATA-SOURCE UPDATE flagged below (stale idr_df.csv -> current per-segment table).
"""IDR density along (A) the protein and (B) the exon, TF vs non-TF.
Interval-level IDR data (idr_df) exists only at canonical/base accession level,
so both panels are per canonical protein (matches the reference figure's per-protein design)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root
NB=100
rel=np.linspace(0,1,NB,endpoint=False)+0.5/NB   # 100 relative-position sample points

# --- IDR intervals per canonical (base) accession ---
# DATA-SOURCE UPDATE: the recovered version read a stale, canonical-only idr_df.csv from an older
# proteome run (28,626 segments, only 33% in the current cohort). This reads the per-segment table
# written by notebooks/00 for the current cohort (>=20 aa); canonical proteins are keyed the same way.
idr=pd.read_csv(REPO/"data/annotated/idr_df.csv", usecols=["isoform_accession","idr_start_aa","idr_end_aa"], low_memory=False)
idr=idr.rename(columns={"isoform_accession":"upkb_accession"})
intervals={}
for acc,s,e in zip(idr["upkb_accession"], idr["idr_start_aa"], idr["idr_end_aa"]):
    intervals.setdefault(str(acc),[]).append((int(s),int(e)))

def profile_from_binary(binary):
    """sample NB relative positions across a 0/1 residue array -> length-NB vector."""
    L=len(binary)
    if L==0: return None
    idx=np.minimum((rel*L).astype(int), L-1)
    return binary[idx].astype(float)

# --- permutation test of terminal (ends vs middle) enrichment, Susie-style ---
# Null: circularly shift each unit's disorder profile (preserves segment structure & amount,
# randomizes position). Statistic = mean(disorder at ends) - mean(disorder in middle).
ENDS=(rel<0.10)|(rel>0.90); MID=(rel>=0.40)&(rel<=0.60)
def _pooled_profile(bins, shifts=None):
    s=np.zeros(NB)
    for i,b in enumerate(bins):
        L=len(b)
        bb=np.roll(b,shifts[i]) if shifts is not None else b
        s+=bb[np.minimum((rel*L).astype(int),L-1)]
    return s/len(bins)
def perm_p(bins, N=200, seed=0):
    from scipy.stats import norm
    rng=np.random.default_rng(seed)
    obs=_pooled_profile(bins); T=obs[ENDS].mean()-obs[MID].mean()
    null=np.empty(N)
    for k in range(N):
        shifts=[int(rng.integers(0,max(1,len(b)))) for b in bins]
        pr=_pooled_profile(bins,shifts); null[k]=pr[ENDS].mean()-pr[MID].mean()
    z=(T-null.mean())/null.std(ddof=1)
    return T, z, float(norm.sf(z))
def fmt_p(z,p):
    import math
    if not np.isfinite(p) or p<1e-30: return "Permutation $P$ < 10$^{-30}$"
    if p>=0.05: return f"Permutation $P$ = {p:.2f} (n.s.)"
    if p<1e-3:
        e=int(math.floor(math.log10(p))); m=p/10**e
        return f"Permutation $P$ = {m:.0f}×10$^{{{e}}}$"
    return f"Permutation $P$ = {p:.3f}"

# ================= PANEL A: along the protein =================
m=pd.read_csv(REPO/"data/annotated/IDRisoforms_df_geneage.csv",
              usecols=["isoform_accession","is_canonical","tf_group","Length"], low_memory=False)
canon=m[m["is_canonical"]==True].copy()
accA={"TF":[], "Non-TF":[]}; binA={"TF":[], "Non-TF":[]}
for acc,grp,L in zip(canon["isoform_accession"], canon["tf_group"], canon["Length"]):
    L=int(L)
    if L<10: continue
    b=np.zeros(L,dtype=np.int8)
    for s,e in intervals.get(str(acc),[]):
        b[max(0,s-1):min(L,e)]=1
    p=profile_from_binary(b)
    if p is not None: accA[str(grp)].append(p); binA[str(grp)].append(b)
profA={g:np.vstack(v).mean(0) for g,v in accA.items()}
nA={g:len(v) for g,v in accA.items()}
TA,zA,pA=perm_p(binA["TF"]+binA["Non-TF"], N=200, seed=0); pAtxt=fmt_p(zA,pA)
print("Panel A (protein) n:", nA, "| perm T=%.3f z=%.1f p=%.1e"%(TA,zA,pA))

# ================= PANEL B: along the exon =================
ex=pd.read_csv(REPO/"data/genomic/processed/exon_level_table.csv",
               usecols=["isoform_accession","is_canonical","is_coding_exon","aa_start","aa_end","is_tf"], low_memory=False)
ex=ex[(ex["is_canonical"]==True) & (ex["is_coding_exon"]==True)].copy()
ex=ex.dropna(subset=["aa_start","aa_end"])
accB={"TF":[], "Non-TF":[]}; binB={"TF":[], "Non-TF":[]}
for acc,a,bb,tf in zip(ex["isoform_accession"], ex["aa_start"], ex["aa_end"], ex["is_tf"]):
    a=int(a); bb=int(bb)
    if bb-a+1 < 6: continue                        # need a few residues
    ivs=intervals.get(str(acc),[])
    Lex=bb-a+1
    b=np.zeros(Lex,dtype=np.int8)
    for s,e in ivs:                                # mark residues of this exon that fall in an IDR
        lo=max(a,s); hi=min(bb,e)
        if hi>=lo: b[lo-a:hi-a+1]=1
    p=profile_from_binary(b)
    if p is not None:
        gB="TF" if bool(tf) else "Non-TF"; accB[gB].append(p); binB[gB].append(b)
profB={g:np.vstack(v).mean(0) for g,v in accB.items()}
nB={g:len(v) for g,v in accB.items()}
allB=binB["TF"]+binB["Non-TF"]; rngB=np.random.default_rng(1)
if len(allB)>6000: allB=[allB[i] for i in rngB.choice(len(allB),6000,replace=False)]
TB,zB,pB=perm_p(allB, N=200, seed=1); pBtxt=fmt_p(zB,pB)
print("Panel B (exon) n:", nB, "| perm T=%.3f z=%.1f p=%.1e"%(TB,zB,pB))

# ================= FIGURE (reference style) =================
x=rel*100
STYLE={"TF":dict(color="#E4761B",ls="-",lw=2.4), "Non-TF":dict(color="#4C72B0",ls="-",lw=2.4)}
fig,axes=plt.subplots(1,2,figsize=(9,3.0))
def panel(ax,prof,n,title,left_lab,right_lab,ptxt):
    for g in ["TF","Non-TF"]:
        base = f"{g}s" if g=="TF" else g
        ax.plot(x,prof[g],label=f"{base} (n = {n[g]:,})",**STYLE[g])
    ax.set_ylim(0,1); ax.set_xlim(0,100)
    ax.set_xticks([0,20,40,60,80,100])
    ax.tick_params(labelsize=8)
    ax.set_title(title,fontsize=11,fontweight="bold")
    ax.text(0.035,0.94,ptxt,transform=ax.transAxes,fontsize=8,va="top")
    tr=ax.get_xaxis_transform()   # x in data coords, y in axes fraction -> labels ON the x-axis
    ax.text(0,-0.18,left_lab,transform=tr,ha="left",va="top",fontsize=8.5)
    ax.text(100,-0.18,right_lab,transform=tr,ha="right",va="top",fontsize=8.5)
panel(axes[0],profA,nA,"IDR Density Along the Protein","N-terminus","C-terminus",pAtxt)
axes[0].set_ylabel("Proportion Disordered",fontsize=10)
axes[0].legend(frameon=False,loc="upper right",fontsize=8.5)
panel(axes[1],profB,nB,"IDR Density Along the Exon","Exon 5′","Exon 3′",pBtxt)
axes[1].legend(frameon=False,loc="upper right",fontsize=8.5)
fig.tight_layout()
fig.savefig(REPO/"figures/proteome/fig_idr_position_protein_exon.pdf",bbox_inches="tight")
print("saved.")
# quick numbers
for lab,prof in [("protein",profA),("exon",profB)]:
    for g in ["TF","Non-TF"]:
        pr=prof[g]
        print(f"  {lab} {g}: ends mean={np.mean([pr[:10].mean(),pr[-10:].mean()]):.2f}  middle(40-60%)={pr[40:60].mean():.2f}")
