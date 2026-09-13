# Recovered from session records (the lab repository's reproducibility/recovered_scripts/remake_fig2.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches.
"""Remake fig2A / fig2B: canonical (left) vs ALTERNATIVE (right), instead of canonical vs all.
Reuses the original split_violin style from notebooks/02_proteome_analysis.ipynb."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, gaussian_kde
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root
outdir = REPO / "figures" / "proteome"

all_isoforms_df = pd.read_csv(REPO / "data/annotated/IDRisoforms_df.csv", low_memory=False)
print("All isoforms:", all_isoforms_df.shape)

def format_p(p):
    if p < 1e-300: return "p < 1e-300"
    elif p < 0.001: return f"p = {p:.2e}"
    else: return f"p = {p:.3f}"

def stars(p):
    if p < 1e-4:  return "****"
    elif p < 1e-3: return "***"
    elif p < 1e-2: return "**"
    elif p < 5e-2: return "*"
    else:          return "ns"

def bracket(ax, x1, x2, y, text, color="black", tick=2.0, lw=1.2, fs=12):
    """Draw a significance bracket from x1..x2 at height y with a centered star/ns label."""
    ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick],
            color=color, lw=lw, clip_on=False, solid_capstyle="butt")
    va = "bottom"
    ax.text((x1 + x2) / 2, y + 0.6, text, ha="center", va=va,
            color=color, fontsize=fs, fontweight="bold", clip_on=False)

def split_violin_canon_vs_alt(plot_df, metric, title, filename=None, title_size=9,
                              ylabel=None, figsize=(8, 5), save=True):
    """Left half = canonical isoforms | Right half = ALTERNATIVE (non-canonical) isoforms."""
    non_tf_canon = plot_df.query("tf_group == 'Non-TF' and is_canonical == True")[metric].dropna()
    tf_canon     = plot_df.query("tf_group == 'TF' and is_canonical == True")[metric].dropna()
    non_tf_alt   = plot_df.query("tf_group == 'Non-TF' and is_canonical == False")[metric].dropna()
    tf_alt       = plot_df.query("tf_group == 'TF' and is_canonical == False")[metric].dropna()

    # Across groups: Non-TF vs TF (within canonical; within alternative)
    _, p_canon_non_tf_vs_tf = mannwhitneyu(non_tf_canon, tf_canon, alternative="two-sided")
    _, p_alt_non_tf_vs_tf   = mannwhitneyu(non_tf_alt, tf_alt, alternative="two-sided")
    # Within groups: canonical vs alternative
    _, p_non_tf_canon_vs_alt = mannwhitneyu(non_tf_canon, non_tf_alt, alternative="two-sided")
    _, p_tf_canon_vs_alt     = mannwhitneyu(tf_canon, tf_alt, alternative="two-sided")

    fig, ax = plt.subplots(figsize=figsize)
    centers = np.array([1, 2])
    xgap = 0.05
    canon_positions = centers - xgap / 2
    alt_positions   = centers + xgap / 2
    violin_width = 0.55
    max_half_width = violin_width / 2

    canon_fill = "#4C72B0"; alt_fill = "#DD8452"
    canon_dark = "#2F4B7C"; alt_dark = "#A95A2C"

    def kde_width_at_y(vals, y, max_half_width):
        vals = np.asarray(vals.dropna())
        if len(vals) < 2 or np.std(vals) == 0:
            return max_half_width * 0.20
        kde = gaussian_kde(vals)
        y_grid = np.linspace(vals.min(), vals.max(), 300)
        max_density = kde(y_grid).max()
        density_at_y = kde([y])[0]
        return max_half_width * (density_at_y / max_density)

    vp_left = ax.violinplot([non_tf_canon, tf_canon], positions=canon_positions,
                            widths=violin_width, showmeans=False, showmedians=False,
                            showextrema=False, side="low")
    for body in vp_left["bodies"]:
        body.set_facecolor(canon_fill); body.set_edgecolor("black")
        body.set_linewidth(1.0); body.set_alpha(0.75)
    vp_right = ax.violinplot([non_tf_alt, tf_alt], positions=alt_positions,
                             widths=violin_width, showmeans=False, showmedians=False,
                             showextrema=False, side="high")
    for body in vp_right["bodies"]:
        body.set_facecolor(alt_fill); body.set_edgecolor("black")
        body.set_linewidth(1.0); body.set_alpha(0.75)

    for x, vals in zip(canon_positions, [non_tf_canon, tf_canon]):
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        for y in [q1, q3]:
            w = kde_width_at_y(vals, y, max_half_width)
            ax.hlines(y, x - w, x, linewidth=1.1, color=canon_dark, linestyles="dashed")
        w = kde_width_at_y(vals, med, max_half_width)
        ax.hlines(med, x - w, x, linewidth=2.4, color=canon_dark)
    for x, vals in zip(alt_positions, [non_tf_alt, tf_alt]):
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        for y in [q1, q3]:
            w = kde_width_at_y(vals, y, max_half_width)
            ax.hlines(y, x, x + w, linewidth=1.1, color=alt_dark, linestyles="dashed")
        w = kde_width_at_y(vals, med, max_half_width)
        ax.hlines(med, x, x + w, linewidth=2.4, color=alt_dark)

    ax.set_xticks(centers)
    ax.set_xticklabels(["Non-TF", "TF"], fontsize=11)
    # per-half sample counts directly under each violin, color-matched to the legend
    for c, n_canon, n_alt in [(centers[0], non_tf_canon, non_tf_alt),
                              (centers[1], tf_canon, tf_alt)]:
        ax.annotate(f"canonical\nn={len(n_canon):,}", xy=(c - 0.16, 0),
                    xycoords=("data", "axes fraction"), xytext=(0, -26),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=8, color=canon_dark, annotation_clip=False)
        ax.annotate(f"alternative\nn={len(n_alt):,}", xy=(c + 0.16, 0),
                    xycoords=("data", "axes fraction"), xytext=(0, -26),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=8, color=alt_dark, annotation_clip=False)
    ax.set_xlim(0.45, 2.55); ax.set_ylim(-2, 134)
    if ylabel: ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=title_size)

    # ---- significance annotations ----
    # Within groups (canonical vs alternative): short bracket over each group center
    bracket(ax, centers[0] - 0.14, centers[0] + 0.14, 104,
            stars(p_non_tf_canon_vs_alt))
    bracket(ax, centers[1] - 0.14, centers[1] + 0.14, 104,
            stars(p_tf_canon_vs_alt))
    # Across groups (Non-TF vs TF): long brackets, colored to match the half being compared
    bracket(ax, canon_positions[0], canon_positions[1], 116,
            stars(p_canon_non_tf_vs_tf), color=canon_dark)
    bracket(ax, alt_positions[0], alt_positions[1], 125,
            stars(p_alt_non_tf_vs_tf), color=alt_dark)

    # ---- legend INSIDE the axes ----
    legend_handles = [Patch(facecolor=canon_fill, edgecolor="black", alpha=0.75, label="Canonical"),
                      Patch(facecolor=alt_fill, edgecolor="black", alpha=0.75, label="Alternative")]
    ax.legend(handles=legend_handles, frameon=True, edgecolor="black", facecolor="white",
              framealpha=1.0, loc="upper left", bbox_to_anchor=(0.015, 0.72),
              borderaxespad=0, fontsize=9)
    # significance key (small, bottom of figure)
    fig.text(0.5, 0.005,
             "ns  p≥0.05      *  p<0.05      **  p<0.01      ***  p<0.001      ****  p<0.0001",
             ha="center", va="bottom", fontsize=7.5, color="#444444")
    fig.subplots_adjust(left=0.11, right=0.97, top=0.88, bottom=0.24)
    if save and filename:
        fig.savefig(outdir / f"{filename}.pdf", bbox_inches="tight")

    print("\n" + "=" * 70); print(title)
    for lbl, v in [("Non-TF canonical", non_tf_canon), ("Non-TF alternative", non_tf_alt),
                   ("TF canonical", tf_canon), ("TF alternative", tf_alt)]:
        print(f"  {lbl:20s} n={len(v):,}  median={v.median():.2f}  mean={v.mean():.2f}")
    print(f"  Across: canon Non-TF vs TF {format_p(p_canon_non_tf_vs_tf)} | alt Non-TF vs TF {format_p(p_alt_non_tf_vs_tf)}")
    print(f"  Within: Non-TF canon vs alt {format_p(p_non_tf_canon_vs_alt)} | TF canon vs alt {format_p(p_tf_canon_vs_alt)}")
    plt.close(fig)

# Fig 2A: all isoforms
split_violin_canon_vs_alt(
    plot_df=all_isoforms_df, metric="pct_idr",
    title="TFs vs. Non-TFs Isoform Disorder Content (canonical vs alternative)",
    filename="fig2A_split_pct_idr_canon_vs_alt", ylabel="%IDR")

# Fig 2B: IDR-containing isoforms only
fig2b_df = all_isoforms_df.query("n_idr_segments > 0").copy()
split_violin_canon_vs_alt(
    plot_df=fig2b_df, metric="pct_idr",
    title="TFs vs. Non-TFs Isoform Disorder Content Among IDR-Containing Isoforms\n(canonical vs alternative)",
    filename="fig2B_split_pct_idr_canon_vs_alt", ylabel="%IDR")

print("\nSaved PDFs to figures/proteome/")
