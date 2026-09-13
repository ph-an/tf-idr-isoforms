"""RECONSTRUCTION of three proteome-level figures whose code was lost (originals: the lab repository's archive/figures_without_code/):

  fig4C_geneage_pctchangeidr.pdf          canonical %IDR by gene-age bin, TF and non-TF (2026-06-30)
  fig5D_tf_families_pctidr_nisoforms.pdf  %IDR per isoform for TF families ranked by isoforms/gene (2026-06-30)
  tf_canonical_vs_isoform_length_delta.pdf TF canonical minus alternative isoform length (2026-06-24)

Every number printed on the archived PDFs is reproduced by these definitions (identified by testing
candidates against them); layout is approximate. Inputs: data/annotated/IDRisoforms_df_geneage.csv.

fig5D family selection, reproduced exactly: TF families other than "Unknown" with >=3 genes, ranked by
median isoforms per gene using pandas' default (unstable) sort; the first 5 and last 2 are shown.
CAUTION: GTF2I-like, IRF, Rel, Paired box and SAND all have median 5, so which three of them appear is
decided by the sort algorithm, not the data.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

REPO = Path(__file__).resolve().parents[2]  # repository root
OUT_MAIN = REPO / "figures" / "proteome"
AGE = ["<100", "100-500", "500-1000", ">1000"]
AGE_LABELS = ["<100 Ma\n(youngest)", "100-500 Ma", "500-1000 Ma", ">1000 Ma\n(oldest)"]
AGE_COLORS = ["#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"]


def fmt_p(p):
    return f"p = {p:.2e}" if p < 0.001 else f"p = {p:.3f}"


def fig4c(df):
    can = df[(df.is_canonical == True) & df.age_bin.notna()]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), sharey=False)
    lines = []
    for ax, grp, col in [(axes[0], "Non-TF", "#2F4B7C"), (axes[1], "TF", "#A95A2C")]:
        d = can[can.tf_group == grp]
        data = [d.loc[d.age_bin == a, "pct_idr"].to_numpy() for a in AGE]
        vp = ax.violinplot(data, positions=range(4), widths=0.8, showextrema=False)
        for body, c in zip(vp["bodies"], AGE_COLORS):
            body.set_facecolor(c); body.set_edgecolor("black"); body.set_alpha(0.85)
        for i, v in enumerate(data):
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.hlines(med, i - 0.2, i + 0.2, color="black", lw=1.6)
            ax.hlines([q1, q3], i - 0.14, i + 0.14, color="black", lw=0.8, ls="--")
            ax.text(i, -6, f"n={len(v):,}", ha="center", fontsize=7, color="0.4")
        ax.set_xticks(range(4), AGE_LABELS, fontsize=8)
        ax.set_ylabel("% IDR (canonical isoform)")
        ax.set_ylim(-9, 103)
        ax.set_title(grp, color=col, fontweight="bold")
        r = spearmanr(d.gene_age, d.pct_idr)
        lines.append(f"{grp:<7} ρ = {r.statistic:.3f},  {fmt_p(r.pvalue)}")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="black") for c in AGE_COLORS]
    leg = axes[1].legend(handles, ["<100 Ma", "100-500 Ma", "500-1000 Ma", ">1000 Ma"], title="Age bin",
                         loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8, title_fontsize=8)
    axes[1].text(1.03, 0.62, "Spearman ρ / p\n" + "\n".join(lines), transform=axes[1].transAxes, fontsize=8, va="top")
    fig.suptitle("Gene age vs. % IDR: opposing trends in TFs and Non-TFs", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_MAIN / "fig4C_geneage_pctchangeidr.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig4C_geneage_pctchangeidr")


def fig5d_alt(df):
    genes = df.drop_duplicates("base_accession")
    tf_genes = genes[(genes.tf_group == "TF") & (genes.tf_family != "Unknown")]
    fams = (tf_genes.groupby("tf_family").agg(genes=("base_accession", "size"), med=("n_isoforms", "median")).reset_index())
    fams = fams[fams.genes >= 3].sort_values("med", ascending=False)   # default quicksort: see module docstring
    top, bottom = fams.head(5), fams.tail(2)
    chosen = top.tf_family.tolist() + bottom.tf_family.tolist()
    nontf = df[df.tf_group == "Non-TF"].pct_idr
    alltf = df[df.tf_group == "TF"].pct_idr
    panels = [(f, df[df.tf_family == f].pct_idr, int(fams.set_index("tf_family").loc[f, "med"])) for f in chosen]
    panels += [("Unknown", df[(df.tf_group == "TF") & (df.tf_family == "Unknown")].pct_idr, None),
               ("Other", df[(df.tf_group == "TF") & ~df.tf_family.isin(chosen + ["Unknown"])].pct_idr, None),
               ("Non-TF", nontf, None)]
    fig, axes = plt.subplots(2, 5, figsize=(15, 5.6))
    bins = np.linspace(0, 100, 21)
    for k, (ax, (name, vals, med)) in enumerate(zip(axes.flat, panels)):
        ref, ref_label = (alltf, "vs all TFs") if name == "Non-TF" else (nontf, "vs Non-TF")
        ax.hist(vals, bins=bins, color="darkred" if name == "Non-TF" else "#4C72B0", alpha=0.85)
        ax.axvline(vals.mean(), color="black", ls="--", lw=1.5)
        ax.axvline(ref.mean(), color="gray", ls=":", lw=1.8)
        p = mannwhitneyu(vals, ref, alternative="two-sided").pvalue
        text = f"n = {len(vals)} isoforms\nmean = {vals.mean():.1f}%\n{ref_label}\n{fmt_p(p)}"
        if med is not None:
            text += f"\nmed iso/gene = {med}"
        ax.text(0.97, 0.95, text, transform=ax.transAxes, ha="right", va="top", fontsize=7)
        ax.set_title(name, fontsize=10)
        if k >= 5:
            ax.set_xlabel("% IDR")
        if k in (0, 5):
            ax.set_ylabel("Number of isoforms")
    axes[0, 0].text(-0.45, 0.5, "Most\nisoforms/gene", transform=axes[0, 0].transAxes, rotation=90, ha="center", va="center", fontsize=8, color="#2F4B7C", fontweight="bold")
    axes[1, 0].text(-0.45, 0.5, "Fewest\nisoforms/gene", transform=axes[1, 0].transAxes, rotation=90, ha="center", va="center", fontsize=8, color="#A95A2C", fontweight="bold")
    fig.suptitle("% Disorder per isoform: TF families ranked by isoform count per gene\n"
                 "Top row = most isoforms/gene  |  Bottom row = fewest isoforms/gene + Unknown / Other / Non-TF", fontsize=11)
    fig.legend([plt.Line2D([0], [0], color="black", ls="--"), plt.Line2D([0], [0], color="gray", ls=":")],
               ["Panel mean", "Non-TF mean (TF panels) / All-TF mean (Non-TF panel)"], loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUT_MAIN / "fig5D_tf_families_pctidr_nisoforms.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig5D_tf_families_pctidr_nisoforms:", chosen)


def length_delta(df):
    tf = df[df.tf_group == "TF"]
    can_len = tf[tf.is_canonical == True].set_index("base_accession").Length
    alt = tf[tf.is_canonical == False].copy()
    alt["delta"] = alt.base_accession.map(can_len) - alt.Length
    lo, hi = np.percentile(alt.delta, [1, 99])
    shown = alt.delta.clip(lo, hi)          # outliers piled into the edge bins, as in the original
    n_clipped = int(((alt.delta < lo) | (alt.delta > hi)).sum())
    med = alt.delta.median()
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.axvspan(lo, 0, color="#F4CCCC", alpha=0.35, lw=0)
    ax.axvspan(0, hi, color="#C9D3EA", alpha=0.35, lw=0)
    ax.hist(shown, bins=60, color="#4C72B0", alpha=0.85, edgecolor="white", lw=0.3)
    ax.axvline(0, color="black", ls="--", lw=1.3)
    ax.axvline(med, color="#D62728", ls="--", lw=1.5, label=f"median = {med:+.0f} aa")
    ax.set_xlim(lo - 30, hi + 30)
    ax.legend(frameon=False, loc="center right", fontsize=10)
    ax.text(0.985, 0.96, f"Canonical longer (Δ > 0): {100 * (alt.delta > 0).mean():.1f}%\n"
                         f"Isoform longer  (Δ < 0): {100 * (alt.delta < 0).mean():.1f}%",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8"))
    ax.text(0.015, 0.97, f"{n_clipped} outliers clipped (p1–p99 shown)", transform=ax.transAxes, ha="left", va="top", fontsize=9, color="0.5")
    ax.set_xlabel("Canonical length − Isoform length (aa)", fontsize=11)
    ax.set_ylabel("Number of isoforms", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"TF isoform lengths relative to their canonical form\n"
                 f"(N = {len(alt):,} alternative isoforms, {tf.base_accession.nunique()} TF genes)", fontsize=12)
    fig.tight_layout()
    out = REPO / "figures" / "proteome" / "tf_canonical_vs_isoform_length_delta.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote tf_canonical_vs_isoform_length_delta")


def main():
    df = pd.read_csv(REPO / "data" / "annotated" / "IDRisoforms_df_geneage.csv", low_memory=False)
    OUT_MAIN.mkdir(parents=True, exist_ok=True)
    fig4c(df)
    fig5d_alt(df)
    length_delta(df)


if __name__ == "__main__":
    main()
