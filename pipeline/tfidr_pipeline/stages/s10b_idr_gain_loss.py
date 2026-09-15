"""
Stage s10b — fate of IDR exons across canonical↔alternative splicing:
CONSERVED vs REMOVED vs ADDED, by TF status.

For every canonical↔alternative isoform pair (the same pairs as the Step-10 splice
events), each IDR exon (idr_frac >= 0.5) is exactly one of:
    conserves_IDR — present in BOTH isoforms (identical coding exon)       (IDR kept)
    removes_IDR   — in the canonical only, spliced out of the alternative  (IDR lost)
    adds_IDR      — in the alternative only                                (IDR gained)
These three PARTITION every IDR-exon involvement per comparison. removed/added
reproduce the Step-10 splice-event counts exactly (10,485 / 2,640).

Unit = per canonical↔alternative comparison (per pair). "Conserved" is inherently
relative to a specific alternative isoform, and the same exon can be conserved vs
one alternative and removed vs another (they overlap at the unique-exon level), so
per-pair is the unit that partitions cleanly. Recomputes from the coordinate map
(reusing classify_splice_events' pair logic) rather than reading splice_events,
because splice_events stores only the differential (removed/added) exons.

Test: 2×3 [TF / non-TF] × [conserved / removed / added] — does the fate
distribution depend on TF status? (χ² + Cramér's V effect size.)

Outputs:
    analysis/idr_gain_loss_by_tf.csv
    figures/genomic/fig15_idr_gain_loss_by_tf.pdf
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from tfidr_pipeline import config as C
from tfidr_pipeline import utils as U
from tfidr_pipeline.mapping.classify_splice_events import classify_pair, idr_overlap_frac

CONS, REM, ADD = "#009E73", "#0072B2", "#E69F00"   # green=conserved, blue=removed, orange=added
MAIN = {"resolved_unique", "resolved_collapsed_equivalent"}
CATS = ["conserves_IDR", "removes_IDR", "adds_IDR"]
_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _fmt_p(p):
    if p <= 1e-300:
        return "P < 10" + "-300".translate(_SUP)
    if p >= 0.01:
        return f"P = {p:.2f}"
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    return f"P = {m:.1f}×10{str(e).translate(_SUP)}"


def load_idr_exon_fate():
    """One row per (canonical↔alt pair, IDR exon), column `category` in CATS.
    Mirrors classify_splice_events; removed/added reproduce splice_events exactly."""
    res = pd.read_parquet(C.ISOFORM_MAPPING / "isoform_transcript_mapping_resolved.parquet")
    res = res[res.resolved_mapping_status.isin(MAIN) & res.ENST_resolved.notna()].copy()
    res["ENST"] = res.ENST_resolved.astype(str).str.split(".").str[0]

    master = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    canon = master.drop_duplicates("isoform_accession").set_index("isoform_accession")
    res["is_canonical"] = res.isoform_accession.map(canon["is_canonical"]).fillna(False)
    istf = master.drop_duplicates("base_accession").set_index("base_accession")["is_tf"]

    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")
    cds = cds.sort_values(["ensembl_transcript_id", "exon_rank"])
    exons, strand_of = {}, {}
    for enst, g in cds.groupby("ensembl_transcript_id"):
        strand_of[enst] = g.strand.iloc[0]
        exons[enst] = [{"iv": (int(s), int(e)), "aa": (int(a0), int(a1))}
                       for s, e, a0, a1 in zip(g.cds_genomic_start, g.cds_genomic_end,
                                               g.aa_start, g.aa_end)]
    idr_by = U.idr_segments_by_isoform()

    rows = []
    genes, alt_isos, n_pairs = set(), set(), 0
    for gene, g in res.groupby("base_accession"):
        cr = g[g.is_canonical]; alts = g[~g.is_canonical]
        if cr.empty or alts.empty:
            continue
        cr = cr.iloc[0]; c_iso, c_enst = cr.isoform_accession, cr.ENST
        if c_enst not in exons:
            continue
        is_tf = bool(istf.get(gene, False))
        genes.add(gene)
        C_ex = exons[c_enst]; strand = strand_of.get(c_enst, "+")
        for a in alts.itertuples():
            if a.ENST not in exons or a.ENST == c_enst:
                continue
            n_pairs += 1; alt_isos.add(a.isoform_accession)
            A_ex = exons[a.ENST]; Aset = {e["iv"] for e in A_ex}
            # removed / added (differential IDR exons) — reproduces splice_events
            for et, iv, aa, present in classify_pair(C_ex, A_ex, strand):
                cat = {"canonical": "removes_IDR", "alternative": "adds_IDR"}.get(present)
                if cat is None:
                    continue
                owner = c_iso if present == "canonical" else a.isoform_accession
                if idr_overlap_frac(aa[0], aa[1], idr_by.get(owner, [])) < 0.5:
                    continue
                rows.append((gene, is_tf, iv[0], iv[1], cat))
            # conserved (IDR exon present in BOTH isoforms)
            for e in C_ex:
                if e["iv"] in Aset:
                    aa = e["aa"]
                    if idr_overlap_frac(aa[0], aa[1], idr_by.get(c_iso, [])) >= 0.5:
                        rows.append((gene, is_tf, e["iv"][0], e["iv"][1], "conserves_IDR"))

    cat = pd.DataFrame(rows, columns=["base_accession", "is_tf", "start", "end", "category"])
    meta = dict(n_genes=len(genes), n_alt_isoforms=len(alt_isos), n_pairs=n_pairs)
    return cat, meta


def summarize(cat):
    rows = []
    for lab, sub in [("ALL", cat), ("TF", cat[cat.is_tf]), ("NON-TF", cat[~cat.is_tf])]:
        vc = sub.category.value_counts()
        con, rem, add = (int(vc.get(k, 0)) for k in CATS)
        n = con + rem + add
        rows.append(dict(cohort=lab, conserved=con, removed=rem, added=add, n=n,
                         pct_conserved=round(con / n * 100, 1),
                         pct_removed=round(rem / n * 100, 1),
                         pct_added=round(add / n * 100, 1),
                         remove_to_add_ratio=round(rem / max(add, 1), 2)))
    return pd.DataFrame(rows)


def test_tf_vs_nontf(cat):
    """2×3: rows TF/non-TF, cols conserved/removed/added — is fate TF-dependent?"""
    t = np.array([[int((cat[cat.is_tf].category == k).sum()) for k in CATS],
                  [int((cat[~cat.is_tf].category == k).sum()) for k in CATS]])
    chi2, p, _, _ = chi2_contingency(t, correction=False)
    n = t.sum(); v = float(np.sqrt(chi2 / (n * (min(t.shape) - 1))))   # Cramér's V
    return t, chi2, p, v


def make_figure(summ, meta, p, v, path_pdf):
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
    fig, ax = plt.subplots(1, 2, figsize=(12.8, 4.6))
    order = ["ALL", "TF", "NON-TF"]
    s = summ.set_index("cohort")

    # Panel a — 100% partition of IDR-exon fate per comparison
    a = ax[0]; a.grid(axis="x", color=GRID, lw=.9); a.set_axisbelow(True)
    y = np.arange(len(order))[::-1]
    for yi, c in zip(y, order):
        pc, pr, pa = s.loc[c, "pct_conserved"], s.loc[c, "pct_removed"], s.loc[c, "pct_added"]
        a.barh(yi, pc, color=CONS, zorder=3)
        a.barh(yi, pr, left=pc, color=REM, zorder=3)
        a.barh(yi, pa, left=pc + pr, color=ADD, zorder=3)
        a.text(pc / 2, yi, f"{pc:.0f}%", ha="center", va="center", color="white",
               fontsize=9.5, fontweight="bold")
        a.text(pc + pr / 2, yi, f"{pr:.0f}%", ha="center", va="center", color="white",
               fontsize=9, fontweight="bold")
        a.text(pc + pr + pa / 2, yi, f"{pa:.0f}%", ha="center", va="center",
               color="white", fontsize=9, fontweight="bold")
    a.set_yticks(y); a.set_yticklabels(order, fontweight="bold")
    a.set_xlim(0, 100); a.set_xlabel("% of IDR exons (per canonical↔alternative comparison)")
    a.legend(handles=[Patch(fc=CONS, label="conserved (in both)"),
                      Patch(fc=REM, label="removed (canonical only)"),
                      Patch(fc=ADD, label="added (alternative only)")],
             frameon=False, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3)
    a.set_title("a", loc="left", fontsize=13, fontweight="bold")

    # Panel b — absolute counts per fate, per cohort (log scale)
    b = ax[1]; b.grid(axis="y", color=GRID, lw=.9); b.set_axisbelow(True)
    x = np.arange(len(order)); w = 0.26
    bars = [(b.bar(x - w, [s.loc[c, "conserved"] for c in order], w, color=CONS, zorder=3, label="conserved")),
            (b.bar(x,     [s.loc[c, "removed"]   for c in order], w, color=REM,  zorder=3, label="removed")),
            (b.bar(x + w, [s.loc[c, "added"]     for c in order], w, color=ADD,  zorder=3, label="added"))]
    for grp in bars:
        for bb in grp:
            b.text(bb.get_x() + bb.get_width() / 2, bb.get_height() * 1.06,
                   f"{int(bb.get_height()):,}", ha="center", fontsize=7.5, color=INK)
    b.set_yscale("log"); b.set_ylim(100, 1e5)
    b.set_xticks(x); b.set_xticklabels(order, fontweight="bold")
    # Keep the axis title attached to panel b so it does not consume the full
    # inter-panel gutter (or appear to sit on top of panel a).
    b.set_ylabel("number of IDR exons (log scale)", labelpad=3)
    b.legend(frameon=False, fontsize=8.5, loc="upper right", ncol=3)
    b.set_title("b", loc="left", fontsize=13, fontweight="bold")

    fig.suptitle("Alternative splicing leaves most IDR exons intact; among those remodeled, "
                 "removal exceeds addition — equally in TFs and non-TFs",
                 x=0.012, ha="left", fontweight="bold", fontsize=12.5)
    ns = " (n.s.)" if p >= 0.05 else ""
    fig.text(0.012, -0.05,
             f"Every IDR exon (idr_frac ≥ 0.5) in a canonical↔alternative comparison is conserved "
             f"(present in both), removed (canonical only), or added (alternative only) — a clean "
             f"partition per comparison. Counts are per comparison (n = {meta['n_genes']:,} genes, "
             f"{meta['n_pairs']:,} canonical↔alternative comparisons); removed/added reproduce the "
             f"Step-10 splice events. b, TF vs non-TF fate distribution tested by χ² on a 2×3 "
             f"(TF/non-TF × conserved/removed/added): {_fmt_p(p)}{ns}, Cramér's V = {v:.03f} "
             f"({'negligible' if v < 0.1 else 'small' if v < 0.3 else 'moderate'} effect).",
             fontsize=7.4, color=MUTED, wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    # tight_layout leaves an unnecessarily wide gutter for panel b's y label.
    # Set the final spacing explicitly after it has measured all text.
    fig.subplots_adjust(wspace=0.12)
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path_pdf.name}")


def main():
    cat, meta = load_idr_exon_fate()
    summ = summarize(cat)
    t, chi2, p, v = test_tf_vs_nontf(cat)

    C.stamp(summ).to_csv(C.ANALYSIS / "idr_gain_loss_by_tf.csv", index=False)
    fig_dir = C.FIGURES; fig_dir.mkdir(parents=True, exist_ok=True)
    make_figure(summ, meta, p, v,
                fig_dir / "fig15_idr_gain_loss_by_tf.pdf")

    print("\n" + "=" * 70)
    print("IDR EXON FATE (conserved / removed / added), per comparison, by TF")
    print("=" * 70)
    print(f"genes={meta['n_genes']:,}  alt-isoforms={meta['n_alt_isoforms']:,}  "
          f"canonical↔alt comparisons={meta['n_pairs']:,}")
    print(summ.to_string(index=False))
    print(f"\nTF vs non-TF fate (2×3):")
    print(f"  TF:     conserved {t[0,0]:,}  removed {t[0,1]:,}  added {t[0,2]:,}")
    print(f"  non-TF: conserved {t[1,0]:,}  removed {t[1,1]:,}  added {t[1,2]:,}")
    print(f"  chi2={chi2:.2f}  P={p:.3f}  Cramér's V={v:.4f}  "
          f"({'n.s.' if p >= 0.05 else 'significant'})")
    print(f"\n[out] analysis/idr_gain_loss_by_tf.csv")


if __name__ == "__main__":
    main()
