"""
Stage s9c — CONTINUOUS exon IDR fraction (not thresholded), alternative vs
constitutive exons, split by TF / non-TF.

Unlike fig14 (which bins exons at idr_frac >= 0.5), this uses the raw per-exon
`idr_frac` — the fraction of the exon's residues predicted disordered — so no
information is lost at the 0.5 cutoff.

Reads the s9 exon table (analysis/exon_idr_annotation.csv); no recompute.

Outputs:
    analysis/exon_idr_fraction_continuous.csv
    figures/genomic/fig16_exon_idr_fraction_continuous.(png|pdf)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, gaussian_kde
from tfidr_pipeline import config as C

FEATURE, BASE = "#E69F00", "#0072B2"     # orange = alternative, blue = constitutive
_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

def _stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "n.s."


def _fmt_p(p):
    if p <= 1e-300:
        return "P < 10" + "-300".translate(_SUP)
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    return f"P = {m:.1f}×10{str(e).translate(_SUP)}"


def summarize(ex):
    rows = []
    for tf_lab, tf_val in [("non-TF", False), ("TF", True)]:
        for usage in ["constitutive", "alternative"]:
            v = ex[(ex.is_tf == tf_val) & (ex.usage == usage)].idr_frac
            rows.append(dict(cohort=tf_lab, usage=usage, n=len(v),
                             mean_idr_frac=round(v.mean(), 4),
                             sem=round(v.std() / np.sqrt(len(v)), 5),
                             median=round(v.median(), 3)))
    return pd.DataFrame(rows)


def main():
    ex = pd.read_csv(C.ANALYSIS / "exon_idr_annotation.csv")
    ex = ex[ex.usage.isin(["alternative", "constitutive"])].copy()
    # deduplicate to UNIQUE physical exons (one per gene+exon), so an exon shared
    # across isoforms is counted once rather than once per isoform.
    n_rows = len(ex)
    ex = (ex.groupby(["gene", "exon_id"])
            .agg(idr_frac=("idr_frac", "mean"), usage=("usage", "first"),
                 is_tf=("is_tf", "first")).reset_index())
    print(f"[dedup] {n_rows:,} exon-rows -> {len(ex):,} unique exons (gene+exon_id)")
    summ = summarize(ex)
    C.stamp(summ).to_csv(C.ANALYSIS / "exon_idr_fraction_continuous.csv", index=False)

    # alt-vs-const Mann-Whitney within each cohort
    pvals = {}
    for tf_lab, tf_val in [("non-TF", False), ("TF", True)]:
        a = ex[(ex.is_tf == tf_val) & (ex.usage == "alternative")].idr_frac
        c = ex[(ex.is_tf == tf_val) & (ex.usage == "constitutive")].idr_frac
        pvals[tf_lab] = mannwhitneyu(a, c, alternative="greater").pvalue

    # ── figure ──
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    INK, MUTED, GRID = "#1a1a1a", "#666", "#E6E6E6"
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "axes.titlesize": 11, "axes.titleweight": "bold",
                         "axes.edgecolor": "#444", "axes.linewidth": 0.9,
                         "text.color": INK, "xtick.color": "#444", "ytick.color": "#444",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2), gridspec_kw=dict(width_ratios=[1, 1.45]))
    s = summ.set_index(["cohort", "usage"])
    cohorts = ["non-TF", "TF"]; x = np.arange(2); w = 0.36

    def add_n_labels(axis):
        """Sample sizes on a second line under each cohort's x-tick label."""
        for xi, c in zip(x, cohorts):
            na = int(s.loc[(c, "alternative"), "n"]); nc = int(s.loc[(c, "constitutive"), "n"])
            axis.text(xi, -0.075, f"alt n={na:,}  ·  const n={nc:,}",
                      transform=axis.get_xaxis_transform(), ha="center", va="top",
                      fontsize=7.5, color=MUTED)

    # Panel a — mean idr_frac (± 95% CI) grouped bars
    a = ax[0]; a.grid(axis="y", color=GRID, lw=.9); a.set_axisbelow(True)
    for k, (usage, col, off) in enumerate([("alternative", FEATURE, -w / 2),
                                           ("constitutive", BASE, w / 2)]):
        vals = [s.loc[(c, usage), "mean_idr_frac"] for c in cohorts]
        errs = [1.96 * s.loc[(c, usage), "sem"] for c in cohorts]
        a.bar(x + off, vals, w, yerr=errs, color=col, zorder=3, capsize=3,
              error_kw=dict(lw=1, ecolor="#333"))
        for xi, v in zip(x + off, vals):
            a.text(xi, v + 0.015, f"{v:.2f}", ha="center", fontsize=9, color=INK)
    for xi, c in zip(x, cohorts):
        top = max(s.loc[(c, "alternative"), "mean_idr_frac"], s.loc[(c, "constitutive"), "mean_idr_frac"])
        a.plot([xi - w / 2, xi + w / 2], [top + 0.07, top + 0.07], color="#333", lw=1)
        a.text(xi, top + 0.085, _stars(pvals[c]), ha="center", fontsize=13, fontweight="bold")
    a.set_xticks(x); a.set_xticklabels(cohorts, fontweight="bold"); add_n_labels(a)
    a.set_ylabel("mean per-exon IDR fraction"); a.set_ylim(0, 1.0)
    a.legend(handles=[Patch(fc=FEATURE, label="alternative exon"),
                      Patch(fc=BASE, label="constitutive exon")],
             frameon=False, fontsize=8.5, loc="upper left")

    # Panel b — split violins (alternative = left half, constitutive = right half).
    # Quartile/median styling mirrors fig4A_geneage_pct_idr: dashed Q1/Q3 and a solid
    # thicker median, each spanning the violin's KDE width at that y, in a darker shade.
    b = ax[1]; b.grid(axis="y", color=GRID, lw=.9); b.set_axisbelow(True)

    XGAP, VIOLIN_W = 0.05, 0.95      # small gap so the halves don't touch
    MAX_HALF = VIOLIN_W / 2

    def _darken(hexc, f=0.45):
        import matplotlib.colors as mc
        r, g, bl = mc.to_rgb(hexc)
        return (r * f, g * f, bl * f)

    def half_violin(arr, pos, side, col):
        vp = b.violinplot([arr], positions=[pos], widths=VIOLIN_W, showmeans=False,
                          showmedians=False, showextrema=False,
                          side="low" if side == "left" else "high")
        for body in vp["bodies"]:
            body.set_facecolor(col); body.set_edgecolor("black")
            body.set_linewidth(1.0); body.set_alpha(.75)
        dark = _darken(col)
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        # width of the violin at each y, so every line stops exactly at the violin edge
        if len(arr) >= 2 and np.std(arr) > 0:
            kde = gaussian_kde(arr)
            maxd = kde(np.linspace(arr.min(), arr.max(), 300)).max()
            wq1, wmed, wq3 = (MAX_HALF * (d / maxd) for d in kde([q1, med, q3]))
        else:
            wq1 = wmed = wq3 = MAX_HALF * 0.20
        for yv, wd in [(q1, wq1), (q3, wq3)]:                       # quartiles: dashed
            lo, hi = (pos - wd, pos) if side == "left" else (pos, pos + wd)
            b.hlines(yv, lo, hi, lw=1.6, color=dark, linestyles=(0, (4, 2)), zorder=5)
        lo, hi = (pos - wmed, pos) if side == "left" else (pos, pos + wmed)
        b.hlines(med, lo, hi, lw=3.2, color=dark, zorder=6)          # median: solid, thicker

    for i, c in enumerate(cohorts):
        sub = ex[ex.is_tf == (c == "TF")]
        half_violin(sub[sub.usage == "alternative"].idr_frac.values, i - XGAP / 2, "left", FEATURE)
        half_violin(sub[sub.usage == "constitutive"].idr_frac.values, i + XGAP / 2, "right", BASE)
    b.set_xticks(x); b.set_xticklabels(cohorts, fontweight="bold"); add_n_labels(b)
    b.set_ylabel("per-exon IDR fraction"); b.set_ylim(-0.05, 1.10)
    # Panel subtitles removed per user request.

    fig.suptitle("Alternative exons have higher IDR fraction than constitutive — in both TFs and non-TFs",
                 x=0.012, ha="left", fontweight="bold", fontsize=13)
    fig.text(0.012, -0.03,
             f"Continuous per-exon idr_frac (fraction of residues predicted disordered), not thresholded; "
             f"one row per UNIQUE exon (deduplicated across isoforms by gene+exon). "
             f"a, mean ± 95% CI; stars = alternative vs constitutive within cohort (Mann-Whitney, one-sided): "
             f"non-TF {_fmt_p(pvals['non-TF'])}, TF {_fmt_p(pvals['TF'])}.  b, split violins (KDE); "
             f"solid line = median, dashed lines = quartiles (Q1, Q3), each spanning the violin's width at that "
             f"value. Distribution is bimodal (exons are mostly fully ordered or fully disordered), so medians "
             f"sit at 0 or ~1 while the means in panel a fall between.",
             fontsize=7.4, color=MUTED, wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    OUT = C.FIGURES; OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig16_exon_idr_fraction_continuous.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig16_exon_idr_fraction_continuous.pdf", bbox_inches="tight")
    plt.close(fig)

    print("=" * 60)
    print("CONTINUOUS exon IDR fraction — alt vs const, by TF")
    print("=" * 60)
    print(summ.to_string(index=False))
    print(f"\nalt vs const (Mann-Whitney, one-sided):")
    for c in cohorts:
        print(f"  {c:<7} {_fmt_p(pvals[c])}  ({_stars(pvals[c])})")
    print(f"\n[fig] fig16_exon_idr_fraction_continuous.png")


if __name__ == "__main__":
    main()