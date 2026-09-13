# Recovered from session records (the lab repository's reproducibility/recovered_scripts/remake_fig4A.py); see docs/PROVENANCE.md.
# Changes from the recovered original: repo-relative paths; no temp-dir previews/caches.
"""Restyle fig4A (% IDR by gene age, Non-TF vs TF split violin) to match the fig2 treatment:
legend INSIDE the axes + significance stars (brackets) over each age bin, instead of the
external p-value text box. Reuses the original split_violin_by_age style."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, gaussian_kde
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root
outdir = REPO / "figures" / "proteome"

df = pd.read_csv(REPO / "data/annotated/IDRisoforms_df_geneage.csv", low_memory=False)

def format_p(p):
    if p < 1e-300: return "p < 1e-300"
    elif p < 0.001: return f"p = {p:.2e}"
    else: return f"p = {p:.3f}"

def stars(p):
    if pd.isna(p):  return "NA"
    if p < 1e-4:  return "****"
    elif p < 1e-3: return "***"
    elif p < 1e-2: return "**"
    elif p < 5e-2: return "*"
    else:          return "ns"

def bracket(ax, x1, x2, y, text, color="black", tick=2.0, lw=1.2, fs=12):
    ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick],
            color=color, lw=lw, clip_on=False, solid_capstyle="butt")
    ax.text((x1 + x2) / 2, y + 0.6, text, ha="center", va="bottom",
            color=color, fontsize=fs, fontweight="bold", clip_on=False)

age_order  = ["<100", "100-500", "500-1000", ">1000"]
age_labels = ["<100\nMa", "100-500\nMa", "500-1000\nMa", ">1000\nMa"]

def split_violin_by_age(plot_df, metric, title, ylabel, filename, figsize=(9, 5.2), save=True):
    plot_df = plot_df.dropna(subset=["age_bin", "tf_group", metric]).copy()
    plot_df["age_bin"] = pd.Categorical(plot_df["age_bin"], categories=age_order, ordered=True)

    fig, ax = plt.subplots(figsize=figsize)
    centers = np.arange(1, len(age_order) + 1)
    xgap = 0.05
    non_tf_positions = centers - xgap / 2
    tf_positions     = centers + xgap / 2
    violin_width = 0.75
    max_half_width = violin_width / 2

    non_tf_fill = "#4C72B0"; tf_fill = "#DD8452"
    non_tf_dark = "#2F4B7C"; tf_dark = "#A95A2C"
    pvals = {}

    def kde_width_at_y(vals, y, mhw):
        vals = np.asarray(pd.Series(vals).dropna())
        if len(vals) < 2 or np.std(vals) == 0: return mhw * 0.20
        kde = gaussian_kde(vals)
        yg = np.linspace(vals.min(), vals.max(), 300)
        return mhw * (kde([y])[0] / kde(yg).max())

    def qmed(x, vals, side, dark):
        vals = pd.Series(vals).dropna()
        if len(vals) == 0: return
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        for y in [q1, q3]:
            w = kde_width_at_y(vals, y, max_half_width)
            ax.hlines(y, x - w if side == "left" else x, x if side == "left" else x + w,
                      linewidth=1.1, color=dark, linestyles="dashed")
        w = kde_width_at_y(vals, med, max_half_width)
        ax.hlines(med, x - w if side == "left" else x, x if side == "left" else x + w,
                  linewidth=2.4, color=dark)

    tick_labels = []
    for i, age in enumerate(age_order):
        nx, tx = non_tf_positions[i], tf_positions[i]
        nv = plot_df.query("age_bin == @age and tf_group == 'Non-TF'")[metric].dropna()
        tv = plot_df.query("age_bin == @age and tf_group == 'TF'")[metric].dropna()
        if len(nv) >= 2:
            for b in ax.violinplot([nv], positions=[nx], widths=violin_width,
                                   showmeans=False, showmedians=False, showextrema=False, side="low")["bodies"]:
                b.set_facecolor(non_tf_fill); b.set_edgecolor("black"); b.set_linewidth(1.0); b.set_alpha(0.75)
            qmed(nx, nv, "left", non_tf_dark)
        if len(tv) >= 2:
            for b in ax.violinplot([tv], positions=[tx], widths=violin_width,
                                   showmeans=False, showmedians=False, showextrema=False, side="high")["bodies"]:
                b.set_facecolor(tf_fill); b.set_edgecolor("black"); b.set_linewidth(1.0); b.set_alpha(0.75)
            qmed(tx, tv, "right", tf_dark)
        pvals[age] = mannwhitneyu(nv, tv, alternative="two-sided")[1] if len(nv) and len(tv) else np.nan
        tick_labels.append(f"{age} Ma\nNon-TF={len(nv):,}\nTF={len(tv):,}")

    ax.set_xticks(centers); ax.set_xticklabels(tick_labels, fontsize=8.5)
    ax.set_xlabel("Gene age bin", labelpad=8); ax.set_ylabel(ylabel)
    ax.set_xlim(0.4, len(age_order) + 0.6); ax.set_ylim(-4, 122)
    ax.set_title(title)

    # significance bracket over each age bin (Non-TF vs TF)
    for i, age in enumerate(age_order):
        bracket(ax, centers[i] - 0.24, centers[i] + 0.24, 104, stars(pvals[age]))

    # legend INSIDE — top-right band (empty above the violins)
    legend_handles = [Patch(facecolor=non_tf_fill, edgecolor="black", alpha=0.75, label="Non-TF"),
                      Patch(facecolor=tf_fill, edgecolor="black", alpha=0.75, label="TF")]
    ax.legend(handles=legend_handles, frameon=True, edgecolor="black", facecolor="white",
              framealpha=1.0, loc="upper right", bbox_to_anchor=(0.988, 0.985),
              borderaxespad=0, fontsize=9, ncol=2, columnspacing=1.0, handlelength=1.3)
    fig.text(0.5, 0.005,
             "ns  p≥0.05      *  p<0.05      **  p<0.01      ***  p<0.001      ****  p<0.0001",
             ha="center", va="bottom", fontsize=7.5, color="#444444")
    fig.subplots_adjust(left=0.09, right=0.985, top=0.90, bottom=0.24)

    if save:
        fig.savefig(outdir / f"{filename}.pdf", bbox_inches="tight")
    print(title)
    for age in age_order:
        print(f"  {age:>9s}: Non-TF vs TF  {format_p(pvals[age])}  -> {stars(pvals[age])}")
    plt.close(fig)

split_violin_by_age(df, "pct_idr",
                    title="% IDR by Gene Age: TFs vs Non-TFs",
                    ylabel="%IDR",
                    filename="fig4A_geneage_pct_idr_stars")
print("\nSaved to figures/proteome/")
