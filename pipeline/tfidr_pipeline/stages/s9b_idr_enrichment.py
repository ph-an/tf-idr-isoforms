"""
Stage s9b — Susie's contingency-table test for IDR enrichment in alternative vs
constitutive exons, run separately for TFs and non-TFs, at TWO levels:

  Level 1 (pooled):     one 2x2 per cohort  -> chi-square test of independence
  Level 2 (gene-level): a 2x2 per gene      -> Fisher's exact, then summarize

2x2 layout:   rows = constitutive / alternative ; cols = non-IDR / IDR
IDR status  = is_idr_exon (idr_frac >= 0.5).
Odds ratio  = (const_nonIDR * alt_IDR) / (const_IDR * alt_nonIDR); OR > 1 means
              ALTERNATIVE exons are more likely to encode an IDR than constitutive.

Consumes the s9 exon table (analysis/exon_idr_annotation.csv) — no recompute.

Outputs (analysis/):
    idr_enrichment_level1_pooled.csv
    idr_enrichment_level2_pergene.csv
    idr_enrichment_level2_summary.csv
Figure (figures/genomic/):
    fig14_idr_enrichment_contingency.(png|pdf)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, binomtest, false_discovery_control
from tfidr_pipeline import config as C

MIN_EXONS = 4   # gene must have >= 4 classifiable coding exons to be testable

OKABE = dict(blue="#0072B2", orange="#E69F00", verm="#D55E00", green="#009E73",
             sky="#56B4E9", gray="#8A8A8A")


# ── shared 2x2 helpers ───────────────────────────────────────────────────────
def load_exons():
    ex = pd.read_csv(C.ANALYSIS / "exon_idr_annotation.csv")
    ex = ex[ex.usage.isin(["alternative", "constitutive"])].copy()
    ex["is_idr_exon"] = ex.is_idr_exon.astype(bool)
    return ex


def table_2x2(sub):
    """rows [constitutive, alternative] x cols [non-IDR, IDR] -> 2x2 count array."""
    const = sub.usage == "constitutive"
    idr = sub.is_idr_exon
    a = int((const & ~idr).sum())   # constitutive, non-IDR
    b = int((const & idr).sum())    # constitutive, IDR
    c = int((~const & ~idr).sum())  # alternative,  non-IDR
    d = int((~const & idr).sum())   # alternative,  IDR
    return np.array([[a, b], [c, d]])


def or_haldane(t):
    """Odds ratio of alt-vs-const being IDR, Haldane-corrected so it is always finite."""
    a, b, c, d = t[0, 0], t[0, 1], t[1, 0], t[1, 1]
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a * d) / (b * c)


# ── Level 1: pooled chi-square per cohort ────────────────────────────────────
def level1(ex):
    rows = []
    for lab, sub in [("all", ex), ("TF", ex[ex.is_tf]), ("non-TF", ex[~ex.is_tf])]:
        t = table_2x2(sub)
        chi2, p, dof, _ = chi2_contingency(t, correction=False)
        n = int(t.sum())
        cramers_v = float(np.sqrt(chi2 / (n * (min(t.shape) - 1))))
        rows.append(dict(
            cohort=lab,
            const_nonIDR=int(t[0, 0]), const_IDR=int(t[0, 1]),
            alt_nonIDR=int(t[1, 0]), alt_IDR=int(t[1, 1]),
            const_IDR_pct=round(t[0, 1] / t[0].sum() * 100, 1),
            alt_IDR_pct=round(t[1, 1] / t[1].sum() * 100, 1),
            chi2=round(chi2, 2), dof=int(dof), p_value=p,
            odds_ratio=round(or_haldane(t), 3),
            cramers_v=round(cramers_v, 4), n_exons=n))
    return pd.DataFrame(rows)


# ── Level 2: per-gene Fisher + summary ───────────────────────────────────────
def level2_pergene(ex, min_exons=MIN_EXONS):
    rows = []
    for gene, g in ex.groupby("gene"):
        if len(g) < min_exons:
            continue
        t = table_2x2(g)
        if (t.sum(axis=1) == 0).any() or (t.sum(axis=0) == 0).any():
            continue  # degenerate table (a row or column margin is zero)
        _, p = fisher_exact(t)                       # exact p, robust to zero cells
        rows.append(dict(
            gene=gene, is_tf=bool(g.is_tf.max()),
            const_nonIDR=int(t[0, 0]), const_IDR=int(t[0, 1]),
            alt_nonIDR=int(t[1, 0]), alt_IDR=int(t[1, 1]),
            n_exons=int(len(g)), odds_ratio=or_haldane(t), fisher_p=p))
    per = pd.DataFrame(rows)
    # difference-in-proportion (bounded, stable effect size for small tables)
    per["const_idr_frac"] = per.const_IDR / (per.const_nonIDR + per.const_IDR)
    per["alt_idr_frac"] = per.alt_IDR / (per.alt_nonIDR + per.alt_IDR)
    per["diff_prop"] = per.alt_idr_frac - per.const_idr_frac
    # BH-FDR within each cohort (TF, non-TF tested separately)
    per["fdr_bh"] = np.nan
    for tf_val in (True, False):
        m = per.is_tf == tf_val
        if m.sum():
            per.loc[m, "fdr_bh"] = false_discovery_control(per.loc[m, "fisher_p"].values, method="bh")
    return per


def level2_summary(per):
    rows = []
    for lab, sub in [("all", per), ("TF", per[per.is_tf]), ("non-TF", per[~per.is_tf])]:
        n = len(sub)
        # direction from raw cross-product (avoids inf/nan on zero cells)
        cross = sub.const_nonIDR * sub.alt_IDR - sub.const_IDR * sub.alt_nonIDR
        n_gt1 = int((cross > 0).sum())   # alt more IDR-enriched
        n_lt1 = int((cross < 0).sum())
        n_dir = n_gt1 + n_lt1
        sign_p = binomtest(n_gt1, n_dir, 0.5, alternative="greater").pvalue if n_dir else np.nan
        rows.append(dict(
            cohort=lab, n_testable=n,
            median_OR=round(float(np.median(sub.odds_ratio)), 3),
            median_diff_prop=round(float(np.median(sub.diff_prop)), 3),
            pct_OR_gt1=round(n_gt1 / n * 100, 1),
            sign_test_p=sign_p,
            pct_fisher_sig=round((sub.fisher_p < 0.05).mean() * 100, 1),
            pct_sig_bh_fdr=round((sub.fdr_bh < 0.05).mean() * 100, 1)))
    return pd.DataFrame(rows)


# ── figure ───────────────────────────────────────────────────────────────────
# consistent 2-colour palette (colourblind-safe): orange = featured group
# (alternative exons / TF), blue = baseline (constitutive exons / non-TF)
FEATURE, BASE = "#E69F00", "#0072B2"
_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "n.s."


def _fmt_p(p):
    if p <= 1e-300:
        return "P < 10" + "-300".translate(_SUP)
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    return f"P = {m:.1f}×10{str(e).translate(_SUP)}"


def _lighten(hexc, f=0.55):
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hexc)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def make_figure(l1, per, l2, path_png, path_pdf):
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
    fig, ax = plt.subplots(1, 3, figsize=(17.0, 5.0))
    order = ["all", "TF", "non-TF"]
    LABELS = {"all": "ALL", "TF": "TF", "non-TF": "NON-TF"}
    l1i = l1.set_index("cohort"); l2i = l2.set_index("cohort")

    # ── Panel a — Level 1: alt vs constitutive IDR% per cohort (pooled χ²) ──
    a = ax[0]; a.grid(axis="y", color=GRID, lw=.9); a.set_axisbelow(True)
    x = np.arange(len(order)); w = 0.36
    for i, c in enumerate(order):
        av, cv = l1i.loc[c, "alt_IDR_pct"], l1i.loc[c, "const_IDR_pct"]
        a.bar(i - w / 2, av, w, color=FEATURE, zorder=3)               # alternative
        a.bar(i + w / 2, cv, w, color=BASE, zorder=3)                  # constitutive
        a.text(i - w / 2, av + 1, f"{av:.0f}", ha="center", fontsize=8.5, color=INK)
        a.text(i + w / 2, cv + 1, f"{cv:.0f}", ha="center", fontsize=8.5, color=INK)
        h = max(av, cv)
        pv, orr, chi = l1i.loc[c, "p_value"], l1i.loc[c, "odds_ratio"], l1i.loc[c, "chi2"]
        a.text(i, h + 15, f"χ² = {chi:,.0f}   OR = {orr:.2f}", ha="center", fontsize=8, color=MUTED)
        a.text(i, h + 9, _stars(pv), ha="center", fontsize=15, color=INK, fontweight="bold")
        a.text(i, h + 5.5, _fmt_p(pv), ha="center", fontsize=7.8, color=MUTED)
    a.set_xticks(x)
    a.set_xticklabels([LABELS[c] for c in order], fontsize=11, fontweight="bold", color=INK)
    a.set_ylabel("% of exons that encode an IDR"); a.set_ylim(0, 90)
    a.legend(handles=[Patch(fc=FEATURE, label="alternative exon"),
                      Patch(fc=BASE, label="constitutive exon")],
             frameon=False, fontsize=8.5, loc="upper left", handlelength=1.1)
    a.set_title("a", loc="left", fontsize=13, fontweight="bold")

    # ── Panel b — Level 2: per-gene Δ IDR fraction (bounded, stable) ──
    b = ax[1]; b.grid(axis="y", color=GRID, lw=.9); b.set_axisbelow(True)
    bins = np.linspace(-1, 1, 33)
    for lab, col in [("non-TF", BASE), ("TF", FEATURE)]:
        d = per[per.is_tf == (lab == "TF")].diff_prop
        sp = l2i.loc[lab, "sign_test_p"]
        b.hist(d, bins=bins, density=True, color=col, alpha=.5, zorder=3,
               label=f"{LABELS[lab]}: {l2i.loc[lab,'pct_OR_gt1']:.0f}% of genes Δ>0  "
                     f"({_stars(sp)}, {_fmt_p(sp)})")
        b.axvline(d.median(), color=col, lw=2.2, ls="--", zorder=4)
    b.axvline(0, color="#333", lw=1.2, zorder=4)
    b.text(-0.02, 0.55, "Δ = 0", transform=b.get_xaxis_transform(), fontsize=8,
           color=MUTED, ha="right", va="center", rotation=90)
    b.set_xlabel("per-gene Δ IDR fraction (alternative − constitutive)")
    b.set_ylabel("density of genes"); b.set_xlim(-1, 1)
    b.legend(frameon=False, fontsize=8.2, loc="upper left")
    b.set_title("b", loc="left", fontsize=13, fontweight="bold")

    # ── Panel c — Level 2: volcano (effect size vs significance) ──
    c = ax[2]; c.grid(color=GRID, lw=.9); c.set_axisbelow(True)
    for lab, col, s, al in [("non-TF", BASE, 9, .30), ("TF", FEATURE, 18, .70)]:
        sub = per[per.is_tf == (lab == "TF")]
        xv = np.log2(sub.odds_ratio.clip(2 ** -4, 2 ** 4))
        yv = -np.log10(sub.fisher_p.clip(lower=1e-25))
        c.scatter(xv, yv, s=s, color=col, alpha=al, edgecolors="none", zorder=3, label=LABELS[lab])
    c.axhline(-np.log10(0.05), color=MUTED, ls="--", lw=1, zorder=2)
    c.text(3.9, -np.log10(0.05) + 0.4, "P = 0.05", fontsize=7.5, color=MUTED, ha="right")
    c.axvline(0, color="#333", lw=1.2, zorder=2)
    c.set_xlabel("per-gene log₂(odds ratio)   —   right → more IDR-enriched")
    c.set_ylabel("−log₁₀(Fisher P)"); c.set_xlim(-4, 4)
    c.legend(frameon=False, fontsize=8.5, loc="upper left", markerscale=1.6, handletextpad=.2)
    c.set_title("c", loc="left", fontsize=13, fontweight="bold")

    fig.suptitle("Intrinsic disorder is enriched in alternatively-spliced exons",
                 x=0.012, ha="left", fontweight="bold", fontsize=14)
    fig.text(0.012, -0.03,
             "a, Pooled χ² test of independence (2×2: exon type × IDR status, df = 1) per cohort; "
             "IDR status = idr_frac ≥ 0.5, OR = (alt IDR odds)/(const IDR odds).   "
             "b, Per-gene difference in IDR fraction (alternative − constitutive) — a bounded, small-sample-stable "
             "effect size; % of genes with Δ>0 and a binomial sign-test.   "
             "c, Per-gene volcano: log₂(OR) vs −log₁₀(Fisher P), one dot per gene (≥4 classifiable exons).   "
             "Significance: * P<0.05, ** P<0.01, *** P<0.001.",
             fontsize=7.2, color=MUTED, wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(path_png, dpi=200, bbox_inches="tight")
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path_png.name}")


def main():
    ex = load_exons()
    print(f"[load] {len(ex):,} classifiable coding exons "
          f"({ex.gene.nunique():,} genes; TF exons {int(ex.is_tf.sum()):,})")

    l1 = level1(ex)
    per = level2_pergene(ex)
    l2 = level2_summary(per)

    C.stamp(l1).to_csv(C.ANALYSIS / "idr_enrichment_level1_pooled.csv", index=False)
    C.stamp(per).to_csv(C.ANALYSIS / "idr_enrichment_level2_pergene.csv", index=False)
    C.stamp(l2).to_csv(C.ANALYSIS / "idr_enrichment_level2_summary.csv", index=False)

    fig_dir = C.FIGURES; fig_dir.mkdir(parents=True, exist_ok=True)
    make_figure(l1, per, l2,
                fig_dir / "fig14_idr_enrichment_contingency.png",
                fig_dir / "fig14_idr_enrichment_contingency.pdf")

    # ── console summary ──
    print("\n" + "=" * 70)
    print("LEVEL 1 — pooled chi-square (rows const/alt x cols nonIDR/IDR)")
    print("=" * 70)
    print(l1[["cohort", "alt_IDR_pct", "const_IDR_pct", "odds_ratio",
              "chi2", "dof", "p_value", "cramers_v", "n_exons"]].to_string(index=False))
    print("\n" + "=" * 70)
    print(f"LEVEL 2 — per-gene Fisher (genes with >= {MIN_EXONS} classifiable exons)")
    print("=" * 70)
    print(l2.to_string(index=False))
    print("\n[out] analysis/idr_enrichment_level1_pooled.csv, _level2_pergene.csv, _level2_summary.csv")


if __name__ == "__main__":
    main()
