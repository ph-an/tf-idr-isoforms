"""RECONSTRUCTION of fig5D_tf_families_pctidr_nisoforms.pdf (2026-06-30), whose code was lost (original: the lab
repository's archive/figures_without_code/): %IDR per isoform for TF families ranked by isoforms per gene.

Every number printed on the archived PDF is reproduced by this definition (identified by testing candidates
against it); layout is approximate. Inputs: data/annotated/IDRisoforms_df_geneage.csv.
The same reconstruction originally also drew fig4C (identical to figures/proteome_extended/fig_age1) and a
TF-only length-delta figure (identical to a panel of fig6C); both were removed as duplicates.

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
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[2]  # repository root
OUT_MAIN = REPO / "figures" / "proteome"


def fmt_p(p):
    return f"p = {p:.2e}" if p < 0.001 else f"p = {p:.3f}"


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


def main():
    df = pd.read_csv(REPO / "data" / "annotated" / "IDRisoforms_df_geneage.csv", low_memory=False)
    OUT_MAIN.mkdir(parents=True, exist_ok=True)
    fig5d_alt(df)


if __name__ == "__main__":
    main()
