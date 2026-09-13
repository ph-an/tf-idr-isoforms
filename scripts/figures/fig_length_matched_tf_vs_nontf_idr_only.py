"""RECONSTRUCTION of figures/exon_density/fig_length_matched_tf_vs_nontf_idr_only.pdf (2026-07-01).

This single panel is panel 3 ("IDR regions only: TF vs Non-TF") of the figure made by cell 31 of
notebooks/06_exon_density_length_matched.ipynb; the code that saved it on its own is not in the notebook.
The statistics below are cell 31's, restated; every n and significance mark printed on the archived
PDF is reproduced.

Exon density = (number of exon chunks - 1) / region bp x 1000 per RefSeq transcript, from the
proteome-wide BED files (chrX/chrY excluded; transcripts need >=100 bp of IDR and of non-IDR).
TF transcripts are those in the Lambert et al. 2018 census (mmc7.xlsx). Bins by IDR bp (right-closed);
two-sided Mann-Whitney U per bin.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[2]  # repository root
SUSIE = REPO / "data" / "raw" / "exon_density"   # Susie Song's RefSeq-based inputs
BED = SUSIE / "bed_all_transcripts"
BINS = [0, 300, 600, 1200, 2400, np.inf]
LABELS = ["<300", "300–600", "600–1,200", "1,200–2,400", ">2,400"]
TF_C, NONTF_C = "#2166AC", "#B2182B"


def read_bed(name):
    d = pd.read_csv(BED / name, sep="\t", header=None, names="chr start end tf strand".split())
    d = d[~d.chr.isin(["chrX", "chrY"])]
    return d.assign(len=d.end - d.start)


def main():
    idr = pd.read_excel(SUSIE / "1-s2.0-S0092867420304815-mmc7.xlsx", sheet_name="IDR_classification")
    tf_nm = set("NM_" + idr["RefseqID_with_IDR_index"].str.split("_", n=2, expand=True).iloc[:, 1])
    i, n = read_bed("all_idr.bed"), read_bed("all_nonidr.bed")
    s = (n.groupby("tf").size().to_frame("nonidr_exons").join(i.groupby("tf").size().to_frame("idr_exons"))
          .join(n.groupby("tf")["len"].sum().to_frame("nonidr_bp").join(i.groupby("tf")["len"].sum().to_frame("idr_bp"))))
    s["idr_norm"] = (s.idr_exons - 1) / s.idr_bp * 1000
    s = s.query("nonidr_bp >= 100 and idr_bp >= 100").copy()
    s["is_tf"] = s.index.isin(tf_nm)
    s["bin"] = pd.cut(s.idr_bp, BINS, labels=LABELS, right=True)

    rows = []
    for b in LABELS:
        tv = s.loc[(s.bin == b) & s.is_tf, "idr_norm"].dropna()
        nv = s.loc[(s.bin == b) & ~s.is_tf, "idr_norm"].dropna()
        p = mannwhitneyu(tv, nv, alternative="two-sided").pvalue if len(tv) >= 3 and len(nv) >= 3 else np.nan
        rows.append((b, len(tv), len(nv), tv.median(), nv.median(), p))
        print(f"{b:<12} n_TF={len(tv):>5} n_nonTF={len(nv):>6} TF med={tv.median():.3f} nonTF med={nv.median():.3f} p={p:.4g}")
    r = pd.DataFrame(rows, columns=["bin", "n_tf", "n_nontf", "med_tf", "med_nontf", "p"])
    star = lambda p: "" if pd.isna(p) else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

    x, w = np.arange(len(LABELS)), 0.35
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.bar(x - w / 2, r.med_tf, w, color=TF_C, alpha=0.85, label="TF IDR regions")
    ax.bar(x + w / 2, r.med_nontf, w, color=NONTF_C, alpha=0.85, label="Non-TF IDR regions")
    for k, row in r.iterrows():
        ax.text(k, max(row.med_tf, row.med_nontf) + 0.15, star(row.p), ha="center", va="bottom", fontsize=12, fontweight="bold")
        ax.text(k - w / 2, -0.55, f"n = {row.n_tf:,}", ha="center", va="top", fontsize=7, color=TF_C, fontweight="bold")
        ax.text(k + w / 2, -0.55, f"n = {row.n_nontf:,}", ha="center", va="top", fontsize=7, color=NONTF_C, fontweight="bold")
    ax.set_xticks(x, LABELS)
    ax.set_xlabel("IDR segment length (bp)", labelpad=18)
    ax.set_ylabel("Median exon density (exons per kb)")
    ax.set_ylim(-1.2, r[["med_tf", "med_nontf"]].to_numpy().max() * 1.04)
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.legend(frameon=True, fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 0.92))
    ax.text(0.98, 0.97, "* p<0.05   ** p<0.01   *** p<0.001", transform=ax.transAxes, ha="right", va="top", fontsize=7, color="gray")
    ax.set_title("IDR regions only: TF vs Non-TF exon density\nwithin matched length bins", fontweight="bold")
    fig.tight_layout()
    out = REPO / "figures" / "exon_density" / "fig_length_matched_tf_vs_nontf_idr_only.pdf"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out.relative_to(REPO))


if __name__ == "__main__":
    main()
