"""
Stage s9d — GENE-LEVEL continuous IDR fraction.

For each gene, compute the mean per-exon idr_frac of its ALTERNATIVE exons and of
its CONSTITUTIVE exons (using unique exons, deduplicated across isoforms). One
paired point per gene. This is the most rigorous form of the alt-vs-constitutive
disorder comparison — it fixes BOTH:
  * the 0.5 threshold  (continuous idr_frac, like fig16), and
  * exon-level pseudoreplication (the gene is the unit, like fig14b).

Because each gene contributes both an alternative mean and a constitutive mean,
the values are PAIRED within gene -> Wilcoxon signed-rank (one-sided), not
Mann-Whitney. Only genes with >=1 alternative AND >=1 constitutive exon qualify.

Reads analysis/exon_idr_annotation.csv; no recompute.

Outputs:
    analysis/gene_idr_fraction_alt_vs_const.csv     (one row per gene)
    figures/genomic/fig17_gene_idr_fraction.(png|pdf)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from tfidr_pipeline import config as C

TF_COL, NONTF_COL = "#E69F00", "#0072B2"     # orange = TF, blue = non-TF
_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "n.s."


def _fmt_p(p):
    if p <= 1e-300:
        return "P < 10" + "-300".translate(_SUP)
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    return f"P = {m:.1f}×10{str(e).translate(_SUP)}"


def build_gene_table():
    ex = pd.read_csv(C.ANALYSIS / "exon_idr_annotation.csv")
    ex = ex[ex.usage.isin(["alternative", "constitutive"])]
    # unique exons first (avoid over-weighting exons shared across isoforms)
    uniq = (ex.groupby(["gene", "exon_id"])
              .agg(idr_frac=("idr_frac", "mean"), usage=("usage", "first"),
                   is_tf=("is_tf", "first")).reset_index())
    # per gene: mean idr_frac of alt exons and of constitutive exons
    g = uniq.groupby(["gene", "usage"]).idr_frac.mean().unstack()
    g = g.dropna(subset=["alternative", "constitutive"])          # need BOTH to pair
    g["is_tf"] = g.index.map(uniq.groupby("gene").is_tf.max())
    g["delta"] = g["alternative"] - g["constitutive"]
    return g.reset_index()


def main():
    g = build_gene_table()
    C.stamp(g).to_csv(C.ANALYSIS / "gene_idr_fraction_alt_vs_const.csv", index=False)

    stats = {}
    for lab, sub in [("all", g), ("TF", g[g.is_tf]), ("non-TF", g[~g.is_tf])]:
        _, p = wilcoxon(sub.alternative, sub.constitutive, alternative="greater")
        stats[lab] = dict(n=len(sub), pct_pos=round((sub.delta > 0).mean() * 100, 1),
                          median_delta=round(sub.delta.median(), 3), p=p)

    # ── figure ──
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, MUTED, GRID = "#1a1a1a", "#666", "#E6E6E6"
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "axes.titlesize": 11, "axes.titleweight": "bold",
                         "axes.edgecolor": "#444", "axes.linewidth": 0.9,
                         "text.color": INK, "xtick.color": "#444", "ytick.color": "#444",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.2))

    # Panel a — paired scatter: constitutive (x) vs alternative (y) mean idr_frac per gene
    a = ax[0]; a.grid(color=GRID, lw=.9); a.set_axisbelow(True)
    for lab, col, s, al in [("non-TF", NONTF_COL, 10, .28), ("TF", TF_COL, 22, .70)]:
        sub = g[g.is_tf == (lab == "TF")]
        a.scatter(sub.constitutive, sub.alternative, s=s, color=col, alpha=al,
                  edgecolors="none", zorder=3, label=f"{lab} (n={len(sub):,})")
    a.plot([0, 1], [0, 1], color="#333", lw=1.3, ls="--", zorder=4)
    a.text(0.97, 0.90, "y = x", fontsize=8.5, color="#333", rotation=45, ha="right")
    a.text(0.03, 0.97, "above line →\nalt more disordered", fontsize=8, color=MUTED, va="top")
    a.set_xlabel("mean IDR fraction — constitutive exons")
    a.set_ylabel("mean IDR fraction — alternative exons")
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.legend(frameon=False, fontsize=8.5, loc="lower right", markerscale=1.4)
    a.set_title("a   Per-gene paired means", loc="left", fontsize=11)

    # Panel b — per-gene Δ distribution, split by TF
    b = ax[1]; b.grid(axis="y", color=GRID, lw=.9); b.set_axisbelow(True)
    bins = np.linspace(-1, 1, 41)
    for lab, col in [("non-TF", NONTF_COL), ("TF", TF_COL)]:
        sub = g[g.is_tf == (lab == "TF")]
        st = stats[lab]
        b.hist(sub.delta, bins=bins, density=True, color=col, alpha=.55, zorder=3,
               label=f"{lab}: {st['pct_pos']:.0f}% of genes Δ>0  "
                     f"({_stars(st['p'])}, {_fmt_p(st['p'])})")
        b.axvline(sub.delta.median(), color=col, lw=2.2, ls="--", zorder=4)
    b.axvline(0, color="#333", lw=1.2, zorder=4)
    b.text(-0.02, 0.6, "Δ = 0", transform=b.get_xaxis_transform(), fontsize=8,
           color=MUTED, ha="right", va="center", rotation=90)
    b.set_xlabel("per-gene Δ mean IDR fraction  (alternative − constitutive)")
    b.set_ylabel("density of genes"); b.set_xlim(-1, 1)
    b.legend(frameon=False, fontsize=8.3, loc="upper left")
    b.set_title("b   Per-gene Δ (Wilcoxon signed-rank, paired)", loc="left", fontsize=11)

    fig.suptitle("Per gene, alternative exons are more disordered than constitutive — in both TFs and non-TFs",
                 x=0.012, ha="left", fontweight="bold", fontsize=13)
    fig.text(0.012, -0.03,
             f"One point per gene (mean per-exon idr_frac of its alternative vs constitutive unique exons); "
             f"genes with ≥1 of each. n = {stats['all']['n']:,} genes ({stats['TF']['n']} TF, {stats['non-TF']['n']:,} non-TF). "
             f"Paired within gene → Wilcoxon signed-rank, one-sided. Significance: * P<0.05, ** P<0.01, *** P<0.001.",
             fontsize=7.4, color=MUTED, wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    OUT = C.FIGURES; OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig17_gene_idr_fraction.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig17_gene_idr_fraction.pdf", bbox_inches="tight")
    plt.close(fig)

    print("=" * 64)
    print("GENE-LEVEL continuous IDR fraction — alt vs const (paired)")
    print("=" * 64)
    for lab in ["all", "TF", "non-TF"]:
        s = stats[lab]
        print(f"  {lab:<7} n={s['n']:>5}  {s['pct_pos']:>5}% genes Δ>0  "
              f"median Δ={s['median_delta']:+.3f}  Wilcoxon {_fmt_p(s['p'])} ({_stars(s['p'])})")
    print(f"\n[fig] fig17_gene_idr_fraction.png")


if __name__ == "__main__":
    main()
