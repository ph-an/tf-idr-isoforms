"""RECONSTRUCTION of the 2026-07-14 preliminary genomic figures fig01, fig02, fig04-fig08, fig10-fig12.

The code that drew these figures was lost (no copy in the repository, notebooks, or surviving session
records). Each definition below was identified by recomputing the numbers printed on the archived
figure (the lab repository's archive/figures_without_code/genomiccoords_preliminary_2026-07-14/) from pipeline tables made
with the IDR definition in force on 2026-07-14 (TFIDR_MIN_IDR_LEN_AA=0); all printed values match.
Run under the current definition (>=20 aa) the IDR-segment-based panels change accordingly.

Sets and definitions
  MAIN isoforms   resolved_mapping_status in {resolved_unique, resolved_collapsed_equivalent} (15,888)
  IDR segments    tfidr_pipeline.utils.load_idr_frame() (honours MIN_IDR_LEN_AA)
  pct_idr, n IDRs isoform master (notebooks' >=20 aa definition, independent of the pipeline setting)
  exons           coding rows of processed/exon_level_table; IDR exon = idr_frac >= 0.5, exons without
                  idr_frac count as non-IDR; frame-symmetric = CDS length divisible by 3
  alternative     exon_usage == "alternative" / "constitutive" (among the gene's mapped transcripts; the
                  paper definition, adopted 2026-09-13). Exons of single-transcript genes are in neither.
                  The archived originals of fig06, fig07 and fig10 used GENCODE-wide `is_constitutive`
                  (alternative = not in every GENCODE coding transcript); those versions are in
                  the lab repository's archive/figures_superseded/gencode_wide_alternative_exons/.
Notes carried from verification: fig01 bars are not nested (15,888 MAIN and 15,889 scaffold-validated
overlap in 15,822); fig10 shows the 12 highest of the families passing the >=50-exon filter.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]  # repository root
sys.path.insert(0, str(REPO / "pipeline"))
from tfidr_pipeline import config as C  # noqa: E402
from tfidr_pipeline import utils as U  # noqa: E402

OUT = C.FIGURES
GREY, SKY, BLUE, GREEN, ORANGE, VERM, PINK = "#8C8C8C", "#56B4E9", "#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7"
MAIN_STATUS = {"resolved_unique", "resolved_collapsed_equivalent"}
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "0.9", "axes.axisbelow": True})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def label_bars(ax, bars, fmt, pad):
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + pad, fmt(b.get_height()), ha="center", va="bottom", fontsize=8)


def load():
    master = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    res = pd.read_parquet(C.ISOFORM_MAPPING / "isoform_transcript_mapping_resolved.parquet")
    scaf = pd.read_parquet(C.INTERIM / "genomic_scaffold_validation.parquet")
    main = set(res.loc[res.resolved_mapping_status.isin(MAIN_STATUS), "isoform_accession"])
    segs = U.load_idr_frame()
    segs = segs[segs.acc.isin(main)].copy()
    ex = pd.read_parquet(C.PROCESSED / "exon_level_table.parquet")
    ex = ex[ex.is_coding_exon == True].copy()
    ex["idr_exon"] = (ex.idr_frac >= 0.5)            # NaN idr_frac -> False
    ex["alternative"] = ex.exon_usage == "alternative"
    ex["constitutive"] = ex.exon_usage == "constitutive"
    ex["symmetric"] = ex.cds_length_bp % 3 == 0
    return master, res, scaf, main, segs, ex


def fig01(master, scaf, main):
    steps = [("Starting\ndataset", len(master), GREY),
             ("≥1 ENST", int((master.n_candidate_ENST.fillna(0) > 0).sum()), SKY),
             ("Exact seq\nmatch", len(scaf), BLUE),
             ("Coordinates\nverified", int(scaf.genomic_scaffold_validated.sum()), GREEN),
             ("MAIN\n(1 transcript)", len(main), ORANGE)]
    fig, ax = plt.subplots(figsize=(6.5, 3.7))
    bars = ax.bar(range(len(steps)), [s[1] for s in steps], color=[s[2] for s in steps], width=0.68)
    total = steps[0][1]
    for b, (_, n, _) in zip(bars, steps):
        ax.text(b.get_x() + b.get_width() / 2, n + total * 0.01, f"{n:,}\n({100 * n / total:.1f}%)", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(steps)), [s[0] for s in steps], fontsize=8.5)
    ax.set_ylabel("isoforms")
    ax.set_ylim(0, total * 1.12)
    ax.grid(axis="x", visible=False)
    ax.set_title("Isoform attrition funnel", loc="left", fontweight="bold")
    save(fig, "fig01_attrition_funnel")


def fig02(master, main, segs):
    mm = master[master.isoform_accession.isin(main)]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.3))
    k = mm.n_idr_segments.clip(upper=6).value_counts().reindex(range(7), fill_value=0)
    axes[0].bar(k.index, k.values, color=BLUE)
    axes[0].set_xticks(range(7), ["0", "1", "2", "3", "4", "5", "6+"])
    axes[0].set(xlabel="IDRs per isoform", ylabel="isoforms", title="IDRs per isoform")
    lengths = segs["IDR len"].clip(upper=400)
    axes[1].hist(lengths, bins=40, color=GREEN)
    med = segs["IDR len"].median()
    axes[1].axvline(med, color=VERM, ls="--")
    axes[1].text(med + 5, axes[1].get_ylim()[1] * 0.9, f"median {med:.0f}aa", color=VERM, fontsize=8)
    axes[1].set(xlabel="IDR length (aa, capped 400)", ylabel="IDR segments", title="IDR length")
    axes[2].hist(mm.pct_idr, bins=40, color=ORANGE)
    mean = mm.pct_idr.mean()
    axes[2].axvline(mean, color=VERM, ls="--")
    axes[2].text(mean + 1.5, axes[2].get_ylim()[1] * 0.9, f"mean {mean:.0f}%", color=VERM, fontsize=8)
    axes[2].set(xlabel="% of protein disordered", ylabel="isoforms", title="%IDR per protein")
    for ax in axes:
        ax.title.set_fontweight("bold"); ax.title.set_fontsize(9); ax.title.set_ha("left"); ax.title.set_x(0)
    fig.suptitle(f"IDR inventory (MAIN set, n={len(mm):,})", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "fig02_idr_inventory")


def fig04(master, segs):
    length = segs.acc.map(master.drop_duplicates("isoform_accession").set_index("isoform_accession").protein_length_aa)
    mid = (segs["IDR start"] + segs["IDR end"]) / 2 / length
    pct = pd.cut(mid, [0, 1 / 3, 2 / 3, 1.0001], labels=["N-third", "middle", "C-third"]).value_counts(normalize=True) * 100
    pct = pct.reindex(["N-third", "middle", "C-third"])
    fig, ax = plt.subplots(figsize=(5, 3.7))
    bars = ax.bar(pct.index, pct.values, color=[BLUE, GREY, ORANGE], width=0.62)
    label_bars(ax, bars, lambda v: f"{v:.1f}%", 0.8)
    ax.set(ylabel="% of IDR segments", ylim=(0, 50))
    ax.grid(axis="x", visible=False)
    ax.set_title("IDR position along the protein (by midpoint)", loc="left", fontweight="bold")
    save(fig, "fig04_idr_position")


def fig05(master, main):
    mm = master[master.isoform_accession.isin(main)]
    g = mm.groupby("tf_group").agg(pct=("pct_idr", "mean"), k=("n_idr_segments", "mean"),
                                  any_=("n_idr_segments", lambda x: 100 * (x > 0).mean())).reindex(["Non-TF", "TF"])
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.3))
    for ax, col, title, fmt, top in [(axes[0], "pct", "mean %IDR", lambda v: f"{v:.1f}%", 70),
                                     (axes[1], "k", "mean IDRs / isoform", lambda v: f"{v:.2f}", 2.4),
                                     (axes[2], "any_", "% with ≥1 IDR", lambda v: f"{v:.0f}%", 120)]:
        bars = ax.bar(["non-TF", "TF"], g[col].values, color=[GREY, VERM], width=0.6)
        label_bars(ax, bars, fmt, top * 0.01)
        ax.set_ylim(0, top)
        ax.grid(axis="x", visible=False)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=9)
    fig.suptitle("Transcription factors are markedly more disordered", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "fig05_tf_vs_nontf")


def fig06(ex):
    rows = []
    for label, d in [("All isoforms", ex), ("TF only", ex[ex.is_tf == True])]:
        rows.append((label, 100 * d.idr_exon[d.alternative].mean(), 100 * d.idr_exon[d.constitutive].mean()))
    x = np.arange(len(rows)); w = 0.38
    fig, ax = plt.subplots(figsize=(5.5, 3.7))
    b1 = ax.bar(x - w / 2, [r[1] for r in rows], w, color=ORANGE, label="alternative exons")
    b2 = ax.bar(x + w / 2, [r[2] for r in rows], w, color=BLUE, label="constitutive exons")
    for bars in (b1, b2):
        label_bars(ax, bars, lambda v: f"{v:.0f}%", 0.8)
    ax.set_xticks(x, [r[0] for r in rows])
    ax.set(ylabel="% of exons that encode an IDR", ylim=(0, 70))
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="x", visible=False)
    ax.set_title("Alternative exons preferentially encode IDRs", loc="left", fontweight="bold")
    save(fig, "fig06_alt_vs_constitutive")


def fig07(ex):
    classes = [("IDR\nexons", ex.idr_exon, BLUE), ("non-IDR\nexons", ~ex.idr_exon, GREY),
               ("alternative", ex.alternative, ORANGE), ("constitutive", ex.constitutive, SKY),
               ("IDR &\nalternative", ex.idr_exon & ex.alternative, GREEN)]
    vals = [100 * ex.symmetric[mask].mean() for _, mask, _ in classes]
    fig, ax = plt.subplots(figsize=(5.8, 3.7))
    bars = ax.bar(range(len(classes)), vals, color=[c[2] for c in classes], width=0.68)
    label_bars(ax, bars, lambda v: f"{v:.1f}%", 0.8)
    ax.set_xticks(range(len(classes)), [c[0] for c in classes])
    ax.set(ylabel="% frame-symmetric (CDS len ÷ 3 = 0)", ylim=(0, 55))
    ax.grid(axis="x", visible=False)
    ax.set_title("Frame symmetry of exon classes", loc="left", fontweight="bold")
    save(fig, "fig07_frame_symmetry")


def fig08(master, segs):
    tf = segs.acc.isin(set(master.loc[master.is_tf == True, "isoform_accession"]))
    props = [("IDR len", "IDR length (aa)", "{:.1f}"), ("kappa", "κ (charge segregation)", "{:.3f}"),
             ("fract_pro", "fraction proline", "{:.3f}"), ("fract_polar", "fraction polar", "{:.3f}"),
             ("FCR", "frac. charged residues", "{:.3f}"), ("fract_aro", "fraction aromatic", "{:.3f}")]
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.4))
    for ax, (col, title, fmt) in zip(axes.flat, props):
        vals = [segs.loc[~tf, col].mean(), segs.loc[tf, col].mean()]
        bars = ax.bar(["non-TF", "TF"], vals, color=[GREY, VERM], width=0.6)
        label_bars(ax, bars, lambda v, f=fmt: f.format(v), max(vals) * 0.01)
        ax.set_ylim(0, max(vals) * 1.25)
        ax.grid(axis="x", visible=False)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=9)
    fig.suptitle("Biophysics of IDRs: TF vs non-TF", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "fig08_idr_biophysics")


def fig10(master, ex):
    fam = master.drop_duplicates("base_accession").set_index("base_accession").tf_family
    alt = ex[ex.alternative & (ex.is_tf == True)]
    g = alt.assign(family=alt.base_accession.map(fam)).groupby("family").idr_exon.agg(["size", "mean"])
    g = g[g["size"] >= 50].sort_values("mean", ascending=False).head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 4.3))
    ax.barh([f"{f}  (n={int(n)})" for f, n in zip(g.index, g["size"])], 100 * g["mean"], color=PINK)
    for i, v in enumerate(100 * g["mean"]):
        ax.text(v + 0.8, i, f"{v:.0f}%", va="center", fontsize=8)
    ax.set(xlabel="% of alternative exons encoding an IDR", xlim=(0, 100))
    ax.grid(axis="y", visible=False)
    ax.set_title("IDR-targeted splicing by TF family (≥50 alt exons)", fontweight="bold")
    save(fig, "fig10_tf_family")


def fig11(master):
    g = (master.groupby(["base_accession", "gene_name"])
               .agg(n=("pct_idr", "size"), lo=("pct_idr", "min"), hi=("pct_idr", "max"), tf=("is_tf", "max"))
               .reset_index())
    g = g[(g.n >= 3) & (g.lo >= 5) & (g.hi <= 95)].assign(rng=lambda d: d.hi - d.lo)
    top = g.sort_values("rng", ascending=False).head(15).iloc[::-1]
    labels = [f"{n} ★" if t else n for n, t in zip(top.gene_name, top.tf)]
    y = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.hlines(y, top.lo, top.hi, color="0.85", lw=4, zorder=1)
    ax.scatter(top.lo, y, color=SKY, s=45, zorder=2, label="min %IDR")
    ax.scatter(top.hi, y, color=VERM, s=45, zorder=2, label="max %IDR")
    ax.set_yticks(y, labels)
    ax.set(xlabel="%IDR across the gene's isoforms", xlim=(0, 100))
    ax.grid(axis="y", visible=False)
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, -0.09), ncol=2, fontsize=8)
    ax.text(0.0, -0.15, "★ = transcription factor", transform=ax.transAxes, ha="left", fontsize=8, color="0.3")
    ax.set_title("Top within-gene IDR rewiring (filtered: ≥3 isoforms, 5–95% IDR)", fontweight="bold")
    save(fig, "fig11_top_idr_rewiring_genes")


def fig12(master, main):
    usage = pd.read_csv(C.ANALYSIS / "splicing_idr_alt_exon_usage.csv").set_index("cohort")
    frac = pd.read_csv(C.ANALYSIS / "splicing_idr_exon_idr_fraction_tf.csv").set_index("cohort")
    sym = pd.read_csv(C.ANALYSIS / "splicing_idr_symmetric_exons.csv")
    sym = sym[sym.scope == "internal_exons"].set_index("cohort")
    mm = master[master.isoform_accession.isin(main)].groupby("tf_group").pct_idr.mean()
    groups, x, w = ["non_TF", "TF"], np.arange(2), 0.38
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))
    a, b, c, d = axes.flat
    b1 = a.bar(x - w / 2, usage.loc[groups, "alt_exon_idr_pct"], w, color=ORANGE, label="alternative exon")
    b2 = a.bar(x + w / 2, usage.loc[groups, "const_exon_idr_pct"], w, color=BLUE, label="constitutive exon")
    for bars in (b1, b2):
        label_bars(a, bars, lambda v: f"{v:.0f}%", 0.8)
    a.set(ylabel="% of exons that encode an IDR", ylim=(0, 75)); a.legend(frameon=False, fontsize=8, loc="upper left")
    a.set_title("A  Alternative exons preferentially encode IDR — in both groups", loc="left", fontweight="bold", fontsize=9)
    vals = frac.loc[groups, "pooled_exon_idr_fraction"]
    bars = b.bar(x, vals, 0.6, color=[GREY, VERM])
    label_bars(b, bars, lambda v: f"{v:.2f}", 0.01)
    b.text(1, vals.iloc[1] + 0.09, f"{vals.iloc[1] / vals.iloc[0]:.2f}×", ha="center", color=VERM, fontweight="bold")
    b.set(ylabel="pooled exon IDR fraction (Σ IDR aa / Σ aa)", ylim=(0, 0.75))
    b.set_title("B  TF exons encode far more disorder per exon", loc="left", fontweight="bold", fontsize=9)
    b1 = c.bar(x - w / 2, sym.loc[groups, "idr_exon_symmetric_pct"], w, color=VERM, label="IDR exon")
    b2 = c.bar(x + w / 2, sym.loc[groups, "ordered_exon_symmetric_pct"], w, color=GREY, label="ordered exon")
    for bars in (b1, b2):
        label_bars(c, bars, lambda v: f"{v:.0f}%", 0.6)
    c.axhline(100 / 3, color="0.5", ls=":", lw=1)
    c.text(1.45, 100 / 3 + 0.6, "random", fontsize=7, color="0.5")
    c.set(ylabel="% frame-symmetric (splice-compatible)", ylim=(0, 52)); c.legend(frameon=False, fontsize=8, loc="upper right")
    c.set_title("C  TF IDR exons are more often clean splice cassettes", loc="left", fontweight="bold", fontsize=9)
    bars = d.bar(x, mm.reindex(["Non-TF", "TF"]).values, 0.6, color=[GREY, VERM])
    label_bars(d, bars, lambda v: f"{v:.0f}%", 0.8)
    d.set(ylabel="mean % of protein disordered", ylim=(0, 65))
    d.set_title("D  TFs are markedly more disordered overall", loc="left", fontweight="bold", fontsize=9)
    for ax in axes.flat:
        ax.set_xticks(x, ["non-TF", "TF"])
        ax.grid(axis="x", visible=False)
    fig.suptitle("Alternative splicing encodes IDRs in TFs vs non-TFs: same mechanism, larger magnitude", fontweight="bold")
    fig.tight_layout()
    save(fig, "fig12_tf_vs_nontf_splicing")


def main():
    master, res, scaf, main_set, segs, ex = load()
    fig01(master, scaf, main_set)
    fig02(master, main_set, segs)
    fig04(master, segs)
    fig05(master, main_set)
    fig06(ex)
    fig07(ex)
    fig08(master, segs)
    fig10(master, ex)
    fig11(master)
    fig12(master, main_set)


if __name__ == "__main__":
    main()
