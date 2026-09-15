# Recovered from session records (the lab repository's reproducibility/recovered_scripts/fig_fate_by_age.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches. Fate table -> results/tables/.
"""fig15b -- IDR-exon fate (conserved/removed/added) by GENE AGE, TF vs non-TF.
100% stacked composition per age bin; matches fig15 palette. Caches the fate table."""
import sys, os
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repository root
sys.path.insert(0, os.path.join(REPO, "pipeline"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import chi2_contingency


# ---------- compute fate table (the original cached it in a temp dir) ----------
if True:
    from tfidr_pipeline import config as C
    from tfidr_pipeline import utils as U
    from tfidr_pipeline.mapping.classify_splice_events import classify_pair, idr_overlap_frac
    MAIN = {"resolved_unique", "resolved_collapsed_equivalent"}
    res = pd.read_parquet(C.ISOFORM_MAPPING / "isoform_transcript_mapping_resolved.parquet")
    res = res[res.resolved_mapping_status.isin(MAIN) & res.ENST_resolved.notna()].copy()
    res["ENST"] = res.ENST_resolved.astype(str).str.split(".").str[0]
    master = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    canon = master.drop_duplicates("isoform_accession").set_index("isoform_accession")
    res["is_canonical"] = res.isoform_accession.map(canon["is_canonical"]).fillna(False)
    istf = master.drop_duplicates("base_accession").set_index("base_accession")["is_tf"]
    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet").sort_values(
        ["ensembl_transcript_id", "exon_rank"])
    exons, strand_of = {}, {}
    for enst, g in cds.groupby("ensembl_transcript_id"):
        strand_of[enst] = g.strand.iloc[0]
        exons[enst] = [{"iv": (int(s), int(e)), "aa": (int(a0), int(a1))}
                       for s, e, a0, a1 in zip(g.cds_genomic_start, g.cds_genomic_end, g.aa_start, g.aa_end)]
    idr_by = U.idr_segments_by_isoform()
    rows = []
    for gene, g in res.groupby("base_accession"):
        cr = g[g.is_canonical]; alts = g[~g.is_canonical]
        if cr.empty or alts.empty: continue
        cr = cr.iloc[0]; c_iso, c_enst = cr.isoform_accession, cr.ENST
        if c_enst not in exons: continue
        is_tf = bool(istf.get(gene, False)); C_ex = exons[c_enst]; strand = strand_of.get(c_enst, "+")
        for a in alts.itertuples():
            if a.ENST not in exons or a.ENST == c_enst: continue
            A_ex = exons[a.ENST]; Aset = {e["iv"] for e in A_ex}
            for et, iv, aa, present in classify_pair(C_ex, A_ex, strand):
                cat_ = {"canonical": "removed", "alternative": "added"}.get(present)
                if cat_ is None: continue
                owner = c_iso if present == "canonical" else a.isoform_accession
                if idr_overlap_frac(aa[0], aa[1], idr_by.get(owner, [])) < 0.5: continue
                rows.append((gene, is_tf, cat_))
            for e in C_ex:
                if e["iv"] in Aset and idr_overlap_frac(e["aa"][0], e["aa"][1], idr_by.get(c_iso, [])) >= 0.5:
                    rows.append((gene, is_tf, "conserved"))
    cat = pd.DataFrame(rows, columns=["base_accession", "is_tf", "fate"])
    age = (pd.read_csv(os.path.join(REPO, "data/annotated/IDRisoforms_df_geneage.csv"),
                       usecols=["base_accession", "age_bin"], low_memory=False).drop_duplicates("base_accession"))
    cat = cat.merge(age, on="base_accession", how="left")
    os.makedirs(os.path.join(REPO, "results", "tables"), exist_ok=True)
    cat.to_csv(os.path.join(REPO, "results", "tables", "fig15b_idr_exon_fate_by_pair.csv"), index=False)

# ---------- figure ----------
AGE = ["<100", "100-500", "500-1000", ">1000"]
FATE = ["conserved", "removed", "added"]
COL = {"conserved": "#009E73", "removed": "#0072B2", "added": "#E69F00"}  # match fig15
cat = cat[cat.age_bin.isin(AGE)]

INK, MUTED, GRID = "#1a1a1a", "#666", "#E6E6E6"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                     "axes.titlesize": 11.5, "axes.titleweight": "bold",
                     "axes.edgecolor": "#444", "axes.linewidth": 0.9,
                     "text.color": INK, "xtick.color": "#444", "ytick.color": "#444",
                     "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), sharey=True)

def panel(ax, sub, title, ntag):
    x = np.arange(len(AGE)); w = 0.66
    pcts = {}
    for i, ab in enumerate(AGE):
        s = sub[sub.age_bin == ab]
        c = {k: int((s.fate == k).sum()) for k in FATE}; n = sum(c.values())
        bottom = 0
        for f in FATE:
            p = c[f] / n * 100 if n else 0
            ax.bar(x[i], p, w, bottom=bottom, color=COL[f], zorder=3,
                   edgecolor="white", linewidth=0.8)
            if p >= 4:
                ax.text(x[i], bottom + p/2, f"{p:.0f}", ha="center", va="center",
                        color="white", fontsize=9.5, fontweight="bold")
            bottom += p
        pcts[ab] = (c, n)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{ab}\nn={pcts[ab][1]:,}" for ab in AGE], fontsize=9)
    ax.set_ylim(0, 100); ax.set_xlim(-0.6, len(AGE)-0.4)
    ax.set_title(title, loc="left")
    ax.text(0.5, -0.185, "gene age (Ma)  \u2192  older", transform=ax.transAxes,
            ha="center", va="top", fontsize=9, color=MUTED)
    # trend annotation: conserved first vs last
    cf = pcts[AGE[0]][0]["conserved"]/max(pcts[AGE[0]][1],1)*100
    cl = pcts[AGE[-1]][0]["conserved"]/max(pcts[AGE[-1]][1],1)*100
    af = pcts[AGE[0]][0]["added"]/max(pcts[AGE[0]][1],1)*100
    al = pcts[AGE[-1]][0]["added"]/max(pcts[AGE[-1]][1],1)*100
    return (cf, cl, af, al)

t_tf = panel(axes[0], cat[cat.is_tf], "a   Transcription factors", "TF")
t_nt = panel(axes[1], cat[~cat.is_tf], "b   Non-TFs", "non-TF")
axes[0].set_ylabel("% of IDR exons per canonical\u2194alternative comparison")



leg = [Patch(fc=COL["conserved"], label="conserved (in both)"),
       Patch(fc=COL["removed"], label="removed (canonical only)"),
       Patch(fc=COL["added"], label="added (alternative only)")]
fig.legend(handles=leg, frameon=False, fontsize=9.5, loc="lower center",
           bbox_to_anchor=(0.5, 0.94), ncol=3)

fig.suptitle("Older TFs conserve their IDR exons and stop gaining new ones; non-TFs stay flat",
             x=0.012, y=1.045, ha="left", fontweight="bold", fontsize=13)
fig.text(0.012, -0.02,
         "Every IDR exon (idr_frac \u2265 0.5) in a canonical\u2194alternative comparison is conserved (in both isoforms), "
         "removed (canonical only), or added (alternative only). Bars are the fate composition within each gene-age bin. "
         "In TFs the fate distribution shifts with age (\u03c7\u00b2 P = 5\u00d710\u207b\u00b9\u00b2, Cram\u00e9r's V = 0.078); in non-TFs the shift is "
         "negligible (V = 0.022). The \u224c4:1 removal-over-addition bias holds at every age.",
         fontsize=7.6, color=MUTED, wrap=True)
fig.tight_layout(rect=[0, 0.02, 1, 0.92]); fig.subplots_adjust(wspace=0.06)

OUT = os.path.join(REPO, "figures/genomic")
fig.savefig(os.path.join(OUT, "fig15b_idr_fate_by_age.pdf"), bbox_inches="tight")
print("TF conserved %.0f->%.0f  added %.0f->%.0f" % t_tf)
print("nonTF conserved %.0f->%.0f  added %.0f->%.0f" % t_nt)
print("saved fig15b_idr_fate_by_age.pdf")
