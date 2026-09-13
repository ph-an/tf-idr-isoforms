"""RECONSTRUCTION of fig13_gene_grouping_alt_fraction and fig13a_altfraction_boxplot.

The code behind fig13 (2026-07-15) is lost, and fig13a was made by cropping fig13's PNG and
repainting its p-value label by hand (see docs/PROVENANCE.md).
The metric below was identified by testing candidate definitions against every statistic
printed on the archived figure; it matches all of them on the 2026-07 pipeline outputs:
TF n=365 / non-TF n=4,165 genes, medians 0.50 / 0.36, means 0.53 / 0.44, MWU p = 1.36e-10.

Metric (per gene): share of rows of data/genomic/analysis/exon_idr_annotation.csv
(one row per mapped isoform x coding exon, genes with >=2 mapped transcripts, written by s9)
whose `usage` is "alternative" rather than "constitutive".
Caveats: rows are isoform-expanded (constitutive exons count once per isoform); with unique
exons the TF shift persists (0.69 vs 0.56, p = 3.7e-11); with GENCODE-wide `is_constitutive`
from exon_level_table there is no TF difference (p = 0.86).

Point subsampling and jitter are display-only and seeded; they cannot match the lost original.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[2]  # repository root
ANALYSIS = REPO / "data" / "genomic" / "analysis"
OUT = REPO / "figures" / "genomic"
TABLES = REPO / "results" / "tables"
TF_C, NT_C = "#D55E00", "#888888"
MAX_POINTS = 400


def p_label(p):
    mantissa, exponent = f"{p:.1e}".split("e")
    return rf"MWU $P$ = {mantissa}$\times$10$^{{{int(exponent)}}}$"


def gene_alt_fraction():
    ex = pd.read_csv(ANALYSIS / "exon_idr_annotation.csv", usecols=["gene", "usage", "is_tf"])
    ex = ex[ex["usage"].isin(["alternative", "constitutive"])]
    return (ex.groupby("gene")
              .agg(alt_frac=("usage", lambda s: (s == "alternative").mean()), is_tf=("is_tf", "max"))
              .reset_index())


def box_panel(ax, nt, tf, p, rng, label):
    for x, vals, col in [(0, nt, NT_C), (1, tf, TF_C)]:
        show = vals if len(vals) <= MAX_POINTS else rng.choice(vals, MAX_POINTS, replace=False)
        ax.scatter(x + rng.uniform(-0.18, 0.18, len(show)), show, s=7, color=col, alpha=0.45, lw=0, zorder=1)
        ax.scatter([x + 0.02], [vals.mean()], marker="D", s=45, color=col, zorder=4)
        ax.text(x + 0.06, vals.mean(), f"mean {vals.mean():.2f}", color=col, va="center", fontsize=9)
    ax.boxplot([nt, tf], positions=[0, 1], widths=0.5, whis=(0, 100), showfliers=False,
               medianprops=dict(color="#222", lw=1.8), boxprops=dict(color="#222", lw=1.8),
               whiskerprops=dict(color="#555"), capprops=dict(color="#555"), zorder=3)
    ax.set_xticks([0, 1], [f"non-TF\n(n={len(nt):,} genes)", f"TF\n(n={len(tf):,} genes)"])
    ax.set_ylabel("alternative-exon fraction (per gene)")
    ax.set_ylim(-0.02, 1.02)
    ax.text(0.5, 1.0, label, ha="center", va="center", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#555"))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    genes = gene_alt_fraction()
    genes.to_csv(TABLES / "fig13_gene_alt_exon_fraction.csv", index=False)
    tf = genes.loc[genes["is_tf"] == True, "alt_frac"].to_numpy()
    nt = genes.loc[genes["is_tf"] != True, "alt_frac"].to_numpy()
    p = mannwhitneyu(tf, nt, alternative="two-sided").pvalue
    print(f"TF n={len(tf)} median={np.median(tf):.2f} mean={tf.mean():.2f} | "
          f"non-TF n={len(nt)} median={np.median(nt):.2f} mean={nt.mean():.2f} | MWU p={p:.2e}")

    # fig13: box panel + ECDF panel
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.4))
    box_panel(a, nt, tf, p, np.random.default_rng(13), f"Mann-Whitney p = {p:.1e}")
    a.set_title("A  Each dot = one gene (grouped by gene)", loc="left", fontsize=10, fontweight="bold")
    for vals, col, lab in [(nt, NT_C, "non-TF"), (tf, TF_C, "TF")]:
        xs = np.sort(vals)
        b.step(xs, np.arange(1, len(xs) + 1) / len(xs), where="post", color=col, lw=2.2,
               label=f"{lab} (median {np.median(vals):.2f})")
        b.scatter([np.median(vals)], [0.5], color=col, s=30, zorder=3)
    b.axhline(0.5, ls=":", color="#555", lw=1)
    b.set_xlim(0, 1)
    b.set_ylim(0, 1)
    b.set_xlabel("alternative-exon fraction (per gene)")
    b.set_ylabel("cumulative fraction of genes")
    b.set_title("B  TF distribution is shifted toward more alt splicing", loc="left", fontsize=10, fontweight="bold")
    b.legend(frameon=False, loc="lower right")
    b.spines[["top", "right"]].set_visible(False)
    b.grid(axis="y", color="0.9")
    fig.suptitle("Grouping by gene: TF genes are more alternatively spliced than non-TF genes",
                 x=0.01, ha="left", fontsize=12.5, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "fig13_gene_grouping_alt_fraction.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "fig13_gene_grouping_alt_fraction.pdf", bbox_inches="tight")
    plt.close(fig)

    # fig13a: panel A alone, with the P label style used on the poster
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    box_panel(ax, nt, tf, p, np.random.default_rng(13), p_label(p))
    fig.tight_layout()
    fig.savefig(OUT / "fig13a_altfraction_boxplot.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
