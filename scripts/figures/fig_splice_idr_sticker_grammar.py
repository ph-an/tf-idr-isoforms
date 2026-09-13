# Recovered from session records (the lab repository's reproducibility/recovered_scripts/aromatic_add_vs_remove.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches;
# constitutive baseline = exon_usage == "constitutive" (mapped transcripts, the paper definition) instead of
# GENCODE-wide is_constitutive (2026-09-13; old version in the lab repository's archive/figures_superseded/).
"""Does alternative splicing preferentially STRIP aromatic (condensate-competent) IDR sequence?
Compare aromatic (F/Y/W) and charged (D/E/K/R) content of IDR exons that are
REMOVED (present in canonical, spliced out) vs ADDED (spliced into an alternative),
with a CONSTITUTIVE IDR-exon baseline. Dedup physical exons; report medians, MWU p, Cliff's delta."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root

AROM = set("FYW"); CHG = set("DEKR")
def arom(s): return sum(c in AROM for c in s)/len(s) if s else np.nan
def chg(s):  return sum(c in CHG  for c in s)/len(s) if s else np.nan

def cliffs_delta(a, b):
    a=np.asarray(a); b=np.asarray(b)
    # P(a>b)-P(a<b) via rank method
    from scipy.stats import rankdata
    n1,n2=len(a),len(b)
    if n1==0 or n2==0: return np.nan
    r=rankdata(np.concatenate([a,b]))
    r1=r[:n1].sum()
    U1=r1-n1*(n1+1)/2
    return 2*U1/(n1*n2)-1

# ---- sequences ----
master = pd.read_csv(REPO/"data/annotated/IDRisoforms_df.csv", low_memory=False)
seq = dict(zip(master["isoform_accession"], master["Sequence"]))
print("isoforms with sequence:", len(seq))

def slice_aa(iso, a, b):
    s = seq.get(iso)
    if not isinstance(s,str): return None
    a=int(a); b=int(b)
    sub = s[max(0,a-1):b]           # 1-based inclusive -> python
    return sub if len(sub)>=5 else None   # need >=5 aa to compute a fraction

# ---- added / removed IDR exons ----
se = pd.read_csv(REPO/"data/genomic/isoform_mapping/splice_events.csv", low_memory=False)
core = se[se["idr_consequence"].isin(["removes_IDR","adds_IDR"])].copy()
core["source_iso"] = np.where(core["idr_consequence"]=="removes_IDR",
                              core["canonical_isoform"], core["alternative_isoform"])
core["cat"] = core["idr_consequence"].map({"removes_IDR":"Removed","adds_IDR":"Added"})
# dedup physical exon sequence
core = core.drop_duplicates(subset=["source_iso","exon_aa_start","exon_aa_end","cat"])
core["pep"] = [slice_aa(i,a,b) for i,a,b in zip(core["source_iso"],core["exon_aa_start"],core["exon_aa_end"])]
core = core[core["pep"].notna()].copy()

# ---- constitutive IDR-exon baseline ----
ex = pd.read_csv(REPO/"data/genomic/processed/exon_level_table.csv", low_memory=False)
con = ex[(ex["exon_usage"]=="constitutive") & (ex["is_idr_exon"]==True)].copy()
con = con.drop_duplicates(subset=["isoform_accession","aa_start","aa_end"])
con["pep"] = [slice_aa(i,a,b) for i,a,b in zip(con["isoform_accession"],con["aa_start"],con["aa_end"])]
con = con[con["pep"].notna()].copy()
con["cat"]="Constitutive"; con["is_tf"]=con["is_tf"].astype(bool)

def build(df, tfcol):
    return pd.DataFrame({
        "cat": df["cat"].values,
        "is_tf": df[tfcol].astype(bool).values,
        "arom": [arom(p) for p in df["pep"]],
        "chg":  [chg(p)  for p in df["pep"]],
        "len":  [len(p)  for p in df["pep"]],
    })
D = pd.concat([build(con,"is_tf"), build(core,"is_tf")], ignore_index=True)
print("\nN per category (deduped physical exons):")
print(D.groupby(["cat","is_tf"]).size().unstack().rename(columns={False:"non-TF",True:"TF"}).to_string())

# ---- stats: Removed vs Added (the money comparison) + vs Constitutive ----
def report(sub, label):
    print(f"\n===== {label} =====")
    for prop in ["arom","chg","len"]:
        med = sub.groupby("cat")[prop].median()
        line = "  ".join(f"{c}={med.get(c,np.nan):.4f}" for c in ["Constitutive","Added","Removed"] if c in med.index)
        print(f"  {prop:5s} median: {line}")
    for prop in ["arom","chg"]:
        r = sub[sub.cat=="Removed"][prop].dropna(); a = sub[sub.cat=="Added"][prop].dropna()
        if len(r) and len(a):
            p = mannwhitneyu(r,a,alternative="two-sided")[1]
            d = cliffs_delta(r,a)
            print(f"  {prop}: Removed vs Added  p={p:.2e}  Cliff's d={d:+.3f}  (n_rem={len(r)}, n_add={len(a)})")

report(D, "ALL")
report(D[~D.is_tf], "non-TF")
report(D[D.is_tf], "TF")
# length control: length-match removed to added by 25aa bins, recompute arom
print("\n--- length-controlled (arom, All, 25aa bins, weighted by matched added n) ---")
sub=D[D.cat.isin(["Removed","Added"])].copy()
sub["lb"]=(sub["len"]//25).clip(upper=8)
for lb,g in sub.groupby("lb"):
    r=g[g.cat=="Removed"]["arom"]; a=g[g.cat=="Added"]["arom"]
    if len(r)>=20 and len(a)>=20:
        p=mannwhitneyu(r,a,alternative="two-sided")[1]
        print(f"  len {int(lb)*25}-{int(lb)*25+24}: Removed med={r.median():.4f} (n={len(r)}) vs Added med={a.median():.4f} (n={len(a)})  p={p:.1e}")

# ---- FIGURE: aromatic (top) + charged (bottom), Constitutive/Added/Removed, split TF/non-TF ----
def stars(p):
    return "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "ns"
order=["Constitutive","Added","Removed"]
colors={"Constitutive":"#9AA0A6","Added":"#4C72B0","Removed":"#DD8452"}
props=[("arom","Aromatic residues\n(% of IDR exon, F/Y/W)"),
       ("chg","Charged residues\n(% of IDR exon, D/E/K/R)")]
fig,axes=plt.subplots(2,2,figsize=(10.5,8.6),sharey="row")
for row,(prop,ylab) in enumerate(props):
    for col,(tf,lab) in enumerate([(False,"Non-TF"),(True,"TF")]):
        ax=axes[row,col]; sub=D[D.is_tf==tf]
        data=[sub[sub.cat==c][prop].dropna().values*100 for c in order]
        bp=ax.boxplot(data,positions=[1,2,3],widths=0.6,patch_artist=True,showfliers=False,
                      medianprops=dict(color="black",linewidth=2))
        for patch,c in zip(bp["boxes"],order):
            patch.set_facecolor(colors[c]); patch.set_alpha(0.85); patch.set_edgecolor("black")
        ax.set_xticks([1,2,3])
        ax.set_xticklabels([f"{c}\nn={len(sub[sub.cat==c]):,}" for c in order],fontsize=8.5)
        if row==0: ax.set_title(lab,fontsize=13,fontweight="bold")
        if col==0: ax.set_ylabel(ylab,fontsize=10.5)
        r=sub[sub.cat=="Removed"][prop].dropna(); a=sub[sub.cat=="Added"][prop].dropna()
        if len(r) and len(a):
            p=mannwhitneyu(r,a,alternative="two-sided")[1]
            y=max(np.percentile(r*100,88),np.percentile(a*100,88))+2
            ax.plot([2,2,3,3],[y,y+0.7,y+0.7,y],color="black",lw=1.2)
            ax.text(2.5,y+1.0,stars(p),ha="center",fontsize=12,fontweight="bold")
        ax.set_ylim(0, ax.get_ylim()[1]*1.05)
fig.suptitle("Splicing-added IDR exons carry aromatic-rich, charge-depleted 'sticker' grammar — strongest in TFs\n"
             "(star = Added vs Removed)", fontsize=12.5, y=0.98)
fig.text(0.5,0.006,"ns p≥0.05   * p<0.05   ** p<0.01   *** p<0.001   **** p<0.0001",
         ha="center",fontsize=7.5,color="#444")
fig.subplots_adjust(bottom=0.11,top=0.90,wspace=0.08,hspace=0.28)
fig.savefig(REPO/"figures/genomic/fig_splice_idr_sticker_grammar.pdf",bbox_inches="tight")
print("\nSaved figure.")
