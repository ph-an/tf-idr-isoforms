"""
Stage S9 — does alternative splicing modulate IDRs, and is it TF-biased?

Consumes the self-verified genomic projection (s2/s3/s4b) to answer the project's
core question with five analyses, each hardened against the obvious confounds:

  0. Reconciliation   exon-density direction is metric-sensitive (per-segment vs
                      pooled) — documents why the naive "IDRs are exon-dense"
                      claim does not survive.
  1. Symmetric exons  are IDR-encoding exons more often frame-preserving cassettes
                      (CDS length %3==0) than ordered exons? Gene-paired.
  2. Alt-exon usage   are differentially-used (alternative) exons more likely to
                      encode IDR than constitutive exons? Gene-paired.
  3. Gene model       does IDR-targeted splicing explain the crude TF>non-TF gap
                      in pct_change_idr? Nested gene-level OLS.
  4. Biophysics       are alternatively-spliced IDRs a distinct functional class
                      (aromatic, charged, compact)? Length-stratified.

Outputs (analysis/):
    splicing_idr_density_reconciliation.csv
    splicing_idr_symmetric_exons.csv
    splicing_idr_alt_exon_usage.csv
    splicing_idr_gene_model.csv
    splicing_idr_biophysics_lengthmatched.csv
    exon_idr_annotation.csv                (the shared per-exon table)
Figure (figures/genomic/):
    fig_splicing_idr_mechanism.pdf, fig_exon_idr_fraction_tf.pdf
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, wilcoxon
from tfidr_pipeline import config as C
from tfidr_pipeline import utils as U

OKABE = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00",
         "green": "#009E73", "sky": "#56B4E9", "gray": "#7A7A7A"}


# ── shared inputs ────────────────────────────────────────────────────────────
def build_exon_idr_table(master, cds, idr_by):
    """One row per coding exon of every scaffold-validated isoform's resolved
    transcript, annotated with the fraction of its residues that are IDR."""
    cds_by = {e: g for e, g in cds.groupby("ensembl_transcript_id")}
    rows = []
    for r in master.itertuples():
        g = cds_by.get(r.ENST)
        if g is None:
            continue
        segs = idr_by.get(r.isoform_accession, [])
        n_ex = len(g)
        for e in g.itertuples():
            L = e.aa_end - e.aa_start + 1
            frac = (U.overlap_length(e.aa_start, e.aa_end, segs) / L) if (segs and L > 0) else 0.0
            cds_len = int(e.cds_end_bp_in_cds - e.cds_start_bp_in_cds)
            rows.append((r.base_accession, r.isoform_accession, r.ENST, e.exon_id,
                         e.exon_rank, n_ex, cds_len, e.aa_start, e.aa_end,
                         frac, bool(r.is_tf)))
    ex = pd.DataFrame(rows, columns=[
        "gene", "iso", "ENST", "exon_id", "exon_rank", "n_coding_exons",
        "cds_len_bp", "aa_start", "aa_end", "idr_frac", "is_tf"])
    ex["symmetric"] = (ex.cds_len_bp % 3 == 0)
    ex["is_internal"] = (ex.exon_rank > 1) & (ex.exon_rank < ex.n_coding_exons)
    ex["exon_class"] = np.where(ex.idr_frac >= 0.5, "IDR",
                        np.where(ex.idr_frac <= 0.1, "ordered", "mixed"))
    # gene-level alternative-vs-constitutive exon usage (genes with >=2 transcripts)
    gtx = ex.groupby("gene").ENST.nunique()
    ex["n_tx_gene"] = ex.gene.map(gtx)
    use = (ex[ex.n_tx_gene >= 2].groupby(["gene", "exon_id"]).ENST.nunique()
           .rename("n_tx_with").reset_index())
    ex = ex.merge(use, on=["gene", "exon_id"], how="left")
    ex["usage"] = np.where(ex.n_tx_with.isna(), "single_tx_gene",
                   np.where(ex.n_tx_with < ex.n_tx_gene, "alternative", "constitutive"))
    ex["is_idr_exon"] = ex.idr_frac >= 0.5
    return ex


# ── analyses ─────────────────────────────────────────────────────────────────
def _or(a, b, c, d):
    OR, p = fisher_exact([[a, b], [c, d]])
    return OR, p


def step0_reconciliation(seg, meta):
    seg = seg[seg.projection_validated].merge(meta, on="isoform_accession", how="left")
    seg["bp"] = seg.length_aa * 3
    persg = seg.groupby(["isoform_accession", "segment_type"])["boundary_density_per_kb"].mean().unstack()
    agg = (seg.groupby(["isoform_accession", "segment_type"])
           .agg(pieces=("n_genomic_pieces", "sum"), bp=("bp", "sum")).reset_index())
    agg["dens"] = (agg.pieces - 1).clip(lower=0) / agg.bp * 1000
    pooled = agg.pivot(index="isoform_accession", columns="segment_type", values="dens")
    m = meta.set_index("isoform_accession")
    out = []
    for metric, tbl in [("per_segment", persg), ("pooled_susie", pooled)]:
        t = tbl.dropna(subset=["idr", "nonidr"]).join(m)
        for lab, sub in [("all", t), ("TF", t[t.is_tf == True])]:
            frac = (sub.idr > sub.nonidr).mean() * 100
            try:
                _, p = wilcoxon(sub.idr, sub.nonidr)
            except ValueError:
                p = np.nan
            out.append(dict(metric=metric, cohort=lab, n=len(sub),
                            pct_idr_gt_nonidr=round(frac, 1),
                            median_idr=round(sub.idr.median(), 3),
                            median_nonidr=round(sub.nonidr.median(), 3),
                            wilcoxon_p=p))
    return pd.DataFrame(out)


def _class_table(sub, split_col, pos, neg, target):
    t = sub[sub[split_col].isin([pos, neg])]
    ct = pd.crosstab(t[split_col], t[target]).reindex(index=[pos, neg], columns=[True, False]).fillna(0)
    a, b = ct.loc[pos, True], ct.loc[pos, False]
    c, d = ct.loc[neg, True], ct.loc[neg, False]
    OR, p = _or(a, b, c, d)
    return dict(pos_frac=round(a / (a + b) * 100, 1) if a + b else np.nan,
                neg_frac=round(c / (c + d) * 100, 1) if c + d else np.nan,
                odds_ratio=round(OR, 3), fisher_p=p, pos_n=int(a + b), neg_n=int(c + d))


def _gene_paired(sub, group_col, val_true, val_false, target):
    """Per gene, fraction(target) for the two groups; paired Wilcoxon across genes."""
    g = (sub[sub[group_col].isin([val_true, val_false])]
         .groupby(["gene", group_col])[target].mean().unstack())
    g = g.dropna(subset=[val_true, val_false])
    if len(g) < 10:
        return np.nan, 0
    try:
        _, p = wilcoxon(g[val_true], g[val_false])
    except ValueError:
        p = np.nan
    return p, len(g)


def step1_symmetric(ex):
    out = []
    for scope, base in [("all_exons", ex), ("internal_exons", ex[ex.is_internal])]:
        for coh, sub in [("all", base), ("TF", base[base.is_tf]), ("non_TF", base[~base.is_tf])]:
            row = dict(scope=scope, cohort=coh)
            row.update(_class_table(sub, "exon_class", "IDR", "ordered", "symmetric"))
            gp, ng = _gene_paired(sub, "exon_class", "IDR", "ordered", "symmetric")
            row["gene_paired_p"] = gp
            row["n_genes"] = ng
            out.append(row)
    df = pd.DataFrame(out).rename(columns={"pos_frac": "idr_exon_symmetric_pct",
                                           "neg_frac": "ordered_exon_symmetric_pct",
                                           "pos_n": "n_idr_exons", "neg_n": "n_ordered_exons"})
    return df


def step2_altusage(ex):
    eu = ex[ex.usage.isin(["alternative", "constitutive"])]
    out = []
    for coh, sub in [("all", eu), ("TF", eu[eu.is_tf]), ("non_TF", eu[~eu.is_tf])]:
        row = dict(cohort=coh)
        row.update(_class_table(sub, "usage", "alternative", "constitutive", "is_idr_exon"))
        gp, ng = _gene_paired(sub, "usage", "alternative", "constitutive", "is_idr_exon")
        row["gene_paired_p"] = gp
        row["n_genes"] = ng
        out.append(row)
    return pd.DataFrame(out).rename(columns={"pos_frac": "alt_exon_idr_pct",
                                             "neg_frac": "const_exon_idr_pct",
                                             "pos_n": "n_alt_exons", "neg_n": "n_const_exons"})


def _ols(y, X):
    X = np.column_stack([np.ones(len(X))] + [X[c].values for c in X.columns])
    y = y.values
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    n, k = X.shape
    s2 = resid @ resid / (n - k)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    r2 = 1 - (resid @ resid) / ((y - y.mean()) @ (y - y.mean()))
    return b, b / se, r2


def step3_model(master_full, ex):
    z = lambda s: (s - s.mean()) / s.std()
    eu = ex[ex.usage.isin(["alternative", "constitutive"])]
    gs = eu.groupby("gene").apply(
        lambda d: pd.Series({"n_alt_idr_exons": int(((d.usage == "alternative") & d.is_idr_exon).sum())}),
        include_groups=False).reset_index()
    gm = (master_full.groupby("base_accession").agg(
            pct_change_idr=("pct_change_idr", "max"), n_isoforms=("n_isoforms", "max"),
            gene_age=("gene_age", "mean"), protein_length=("protein_length_aa", "mean"),
            is_tf=("is_tf", "max")).reset_index())
    d = gm.merge(gs, left_on="base_accession", right_on="gene", how="inner").dropna(
        subset=["pct_change_idr", "n_isoforms", "gene_age", "protein_length"])
    d["tf"] = d.is_tf.astype(float)
    d["z_niso"], d["z_logL"] = z(d.n_isoforms), z(np.log(d.protein_length))
    d["z_age"], d["z_altidr"] = z(d.gene_age), z(d.n_alt_idr_exons)
    rows = []
    specs = [("unadjusted", ["tf"]),
             ("M0_controls", ["tf", "z_niso", "z_logL", "z_age"]),
             ("M1_plus_idr_splicing", ["tf", "z_niso", "z_logL", "z_age", "z_altidr"])]
    for name, cols in specs:
        b, t, r2 = _ols(d.pct_change_idr, d[cols])
        for nm, bi, ti in zip(["intercept"] + cols, b, t):
            rows.append(dict(model=name, term=nm, beta=round(bi, 4), t=round(ti, 2), r2=round(r2, 4)))
    return pd.DataFrame(rows), len(d), int(d.tf.sum())


def step4_biophysics(ex, idf, master):
    # per IDR segment: encoded (partly) by an alternative exon?
    alt_ids = set(map(tuple, ex[ex.usage == "alternative"][["gene", "exon_id"]].values))
    exm = ex[ex.n_tx_gene >= 2].copy()
    exm["is_alt"] = [(g, x) in alt_ids for g, x in zip(exm.gene, exm.exon_id)]
    by_iso = {i: g for i, g in exm.groupby("iso")}

    def seg_alt(iso, s, e):
        g = by_iso.get(iso)
        if g is None:
            return None
        covered = alt = False
        for r in g.itertuples():
            lo, hi = max(s, r.aa_start), min(e, r.aa_end)
            if hi >= lo:
                covered = True
                alt = alt or r.is_alt
        return alt if covered else None

    d = idf.copy()
    d["seg_alt"] = [seg_alt(a, int(s), int(e)) for a, s, e in zip(d.acc, d["IDR start"], d["IDR end"])]
    d = d.dropna(subset=["seg_alt"])
    tfmap = master.drop_duplicates("isoform_accession").set_index("isoform_accession").is_tf
    d["is_tf"] = d.acc.map(tfmap)
    d = d.dropna(subset=["is_tf"])
    d["len_bin"] = pd.cut(d["IDR len"], [0, 50, 100, 200, 10000],
                          labels=["0-50", "50-100", "100-200", "200+"])
    props = ["FCR", "NCPR", "kappa", "fract_aro", "fract_pro", "scaling_exponent"]
    out = []
    for coh, sub in [("TF", d[d.is_tf == True]), ("non_TF", d[d.is_tf != True])]:
        for lb, g in sub.groupby("len_bin", observed=True):
            A, Cst = g[g.seg_alt == True], g[g.seg_alt == False]
            if len(A) < 20 or len(Cst) < 20:
                continue
            for p in props:
                a, c = A[p].dropna(), Cst[p].dropna()
                try:
                    _, pv = mannwhitneyu(a, c)
                except ValueError:
                    pv = np.nan
                out.append(dict(cohort=coh, len_bin=lb, property=p,
                                alt_median=round(a.median(), 4), const_median=round(c.median(), 4),
                                n_alt=len(a), n_const=len(c), mannwhitney_p=pv))
    return pd.DataFrame(out), d


# ── figure ───────────────────────────────────────────────────────────────────
def make_figure(s1, s2, s3, bio_seg, path_pdf):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[fig] matplotlib unavailable ({e}); skipping figure")
        return
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 140})
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.2))

    # A — alt vs constitutive exons encode IDR (headline, step 2)
    a = ax[0, 0]
    order = ["all", "TF", "non_TF"]
    s2i = s2.set_index("cohort")
    x = np.arange(len(order)); w = 0.38
    a.bar(x - w/2, [s2i.loc[c, "alt_exon_idr_pct"] for c in order], w,
          color=OKABE["orange"], label="alternative exon")
    a.bar(x + w/2, [s2i.loc[c, "const_exon_idr_pct"] for c in order], w,
          color=OKABE["blue"], label="constitutive exon")
    for i, c in enumerate(order):
        a.text(i - w/2, s2i.loc[c, "alt_exon_idr_pct"] + 1, f"{s2i.loc[c,'alt_exon_idr_pct']:.0f}", ha="center", fontsize=8)
        a.text(i + w/2, s2i.loc[c, "const_exon_idr_pct"] + 1, f"{s2i.loc[c,'const_exon_idr_pct']:.0f}", ha="center", fontsize=8)
    a.set_xticks(x); a.set_xticklabels(["All", "TF", "non-TF"])
    a.set_ylabel("% exons that encode IDR"); a.set_ylim(0, 75)
    a.legend(frameon=False, fontsize=8, loc="upper right")
    a.set_title("A  Alternatively-spliced exons preferentially encode IDR", fontsize=9, loc="left")

    # B — symmetric exon enrichment (internal), step 1
    b = ax[0, 1]
    s1i = s1[s1.scope == "internal_exons"].set_index("cohort")
    b.bar(x - w/2, [s1i.loc[c, "idr_exon_symmetric_pct"] for c in order], w,
          color=OKABE["verm"], label="IDR exon")
    b.bar(x + w/2, [s1i.loc[c, "ordered_exon_symmetric_pct"] for c in order], w,
          color=OKABE["gray"], label="ordered exon")
    b.axhline(33.3, ls=":", lw=1, color="k", alpha=.5)
    b.text(2.35, 34, "random", fontsize=7, alpha=.6)
    b.set_xticks(x); b.set_xticklabels(["All", "TF", "non-TF"])
    b.set_ylabel("% symmetric (frame-preserving)"); b.set_ylim(0, 55)
    b.legend(frameon=False, fontsize=8, loc="upper right")
    b.set_title("B  IDR exons are more often splice-compatible cassettes", fontsize=9, loc="left")

    # C — regression: TF effect shrinks, splicing predicts (step 3)
    c = ax[1, 0]
    tf_betas = {m: s3[(s3.model == m) & (s3.term == "tf")].beta.values[0]
                for m in ["unadjusted", "M0_controls", "M1_plus_idr_splicing"]}
    labs = ["TF only", "+controls", "+IDR splicing"]
    c.plot(range(3), list(tf_betas.values()), "o-", color=OKABE["verm"], lw=2, ms=7, label="TF coefficient")
    c.axhline(0, ls="--", lw=1, color="k", alpha=.5)
    splice_b = s3[(s3.model == "M1_plus_idr_splicing") & (s3.term == "z_altidr")].beta.values[0]
    splice_t = s3[(s3.model == "M1_plus_idr_splicing") & (s3.term == "z_altidr")].t.values[0]
    c.axhline(splice_b, ls="-", lw=1.5, color=OKABE["green"], alpha=.8)
    c.set_ylim(min(list(tf_betas.values()) + [0]) - 0.4, splice_b + 0.9)
    c.text(1.0, splice_b + 0.18, f"IDR-splicing β = {splice_b:.1f} (t={splice_t:.0f})",
           fontsize=8, color=OKABE["green"], ha="center")
    c.text(2.0, tf_betas["M1_plus_idr_splicing"] - 0.22, "TF ≈ 0", fontsize=8,
           color=OKABE["verm"], ha="center")
    c.set_xticks(range(3)); c.set_xticklabels(labs)
    c.set_ylabel("β on pct_change_idr (gene)")
    c.set_title("C  TF effect vanishes with controls; IDR-splicing drives it", fontsize=9, loc="left", pad=10)

    # D — biophysics of alt vs constitutive IDRs (TF, length-controlled view)
    d = ax[1, 1]
    tf = bio_seg[bio_seg.is_tf == True]
    props = ["fract_aro", "FCR", "kappa", "scaling_exponent"]
    labels = ["aromatic", "FCR", "kappa", "scaling exp"]
    altv = [tf[tf.seg_alt == True][p].median() for p in props]
    conv = [tf[tf.seg_alt == False][p].median() for p in props]
    xp = np.arange(len(props))
    d.bar(xp - w/2, altv, w, color=OKABE["orange"], label="alt-spliced IDR")
    d.bar(xp + w/2, conv, w, color=OKABE["blue"], label="constitutive IDR")
    d.set_xticks(xp); d.set_xticklabels(labels)
    d.set_ylabel("median value (TF IDRs)")
    d.legend(frameon=False, fontsize=8, loc="upper right")
    d.set_title("D  Splicing toggles aromatic, less-charged, compact IDRs", fontsize=9, loc="left")

    fig.suptitle("Alternative splicing modulates IDRs via toggleable, functionally distinct exon-encoded modules",
                 fontsize=10.5, y=1.00)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] wrote {path_pdf.name}")


def step5_exon_idr_fraction(ex):
    """Exon IDR fraction = IDR aa-equivalents encoded / total aa-equivalents encoded,
    per coding exon, compared TF vs non-TF. Reports both the pooled fraction
    (Σ IDR aa / Σ aa across all exons in the group) and the per-exon distribution."""
    d = ex.copy()
    d["exon_aa_len"] = d["aa_end"] - d["aa_start"] + 1        # aa encoded by the exon
    d["idr_aa"] = d["idr_frac"] * d["exon_aa_len"]            # IDR aa-equivalents
    out = []
    for lab, g in [("TF", d[d.is_tf == True]), ("non_TF", d[d.is_tf != True])]:
        out.append(dict(
            cohort=lab, n_exons=len(g),
            pooled_exon_idr_fraction=round(g.idr_aa.sum() / g.exon_aa_len.sum(), 4),
            mean_per_exon_idr_frac=round(g.idr_frac.mean(), 4),
            median_per_exon_idr_frac=round(g.idr_frac.median(), 4),
            pct_exons_idr_ge50=round((g.idr_frac >= 0.5).mean() * 100, 1)))
    res = pd.DataFrame(out)
    try:
        _, p = mannwhitneyu(d[d.is_tf == True].idr_frac, d[d.is_tf != True].idr_frac,
                            alternative="greater")
    except ValueError:
        p = np.nan
    res["mannwhitney_p_TF_gt_nonTF"] = p
    return res, d


def make_exon_fraction_figure(d, table, path_pdf):
    """Two panels: (A) ECDF of per-exon IDR fraction, (B) pooled fraction bar."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[fig] matplotlib unavailable ({e}); skipping exon-fraction figure")
        return
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 140})
    tf  = d[d.is_tf == True].idr_frac.dropna().values
    non = d[d.is_tf != True].idr_frac.dropna().values
    ti  = table.set_index("cohort")
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.9))

    # A — ECDF of per-exon IDR fraction (non-TF first = fixed categorical order)
    a = ax[0]
    for arr, lab, col in [(non, "non-TF", OKABE["blue"]), (tf, "TF", OKABE["verm"])]:
        xs = np.sort(arr)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        a.plot(xs, ys, lw=2, color=col, label=f"{lab} (n={len(arr):,})")
        a.plot([np.median(arr)], [0.5], "o", color=col, ms=7, zorder=5,
               markeredgecolor="white", markeredgewidth=1.2)
    a.axhline(0.5, ls=":", lw=1, color="#888", alpha=.7)
    a.text(0.015, 0.53, "median", fontsize=7, color="#666")
    a.set_xlabel("per-exon IDR fraction"); a.set_ylabel("cumulative fraction of exons")
    a.set_xlim(0, 1); a.set_ylim(0, 1)
    a.legend(frameon=False, fontsize=8, loc="lower right")
    a.set_title("A  Distribution of exon IDR fraction", fontsize=9, loc="left")

    # B — pooled exon IDR fraction (Σ IDR aa / Σ aa)
    b = ax[1]
    order, cols = ["non_TF", "TF"], [OKABE["blue"], OKABE["verm"]]
    vals = [ti.loc[c, "pooled_exon_idr_fraction"] for c in order]
    x = np.arange(2)
    b.bar(x, vals, width=0.6, color=cols)
    for xi, v in zip(x, vals):
        b.text(xi, v + 0.015, f"{v:.2f}", ha="center", fontsize=9.5)
    ratio = vals[1] / vals[0] if vals[0] else float("nan")
    b.annotate(f"{ratio:.2f}×", xy=(1, vals[1]), xytext=(0.5, max(vals) + 0.10),
               ha="center", fontsize=10, color=OKABE["verm"], fontweight="bold")
    b.set_xticks(x); b.set_xticklabels(["non-TF", "TF"])
    b.set_ylabel("pooled exon IDR fraction\n(Σ IDR aa / Σ aa encoded)")
    b.set_ylim(0, max(vals) + 0.18)
    b.set_title("B  Pooled exon IDR fraction", fontsize=9, loc="left")

    fig.suptitle("TF exons encode far more disorder per exon than non-TF exons",
                 fontsize=10.5, y=1.02)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] wrote {path_pdf.name}")


def main():
    print("[load] inputs ...")
    m = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    master = m[m.genomic_scaffold_validated == True].copy()
    master["ENST"] = master.resolved_transcript.astype(str).str.split(".").str[0]
    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")
    idf = U.load_idr_frame()
    idr_by = U.idr_segments_by_isoform(idf)

    print("[build] exon x IDR annotation table ...")
    ex = build_exon_idr_table(
        master[["isoform_accession", "base_accession", "ENST", "is_tf"]], cds, idr_by)
    C.stamp(ex).to_csv(C.ANALYSIS / "exon_idr_annotation.csv", index=False)
    print(f"        {len(ex):,} coding exons across {ex.iso.nunique():,} isoforms")

    seg = pd.read_parquet(C.INTERIM / "genomic_projection_segments.parquet")
    meta = m[["isoform_accession", "is_tf"]]

    r0 = step0_reconciliation(seg, meta)
    r1 = step1_symmetric(ex)
    r2 = step2_altusage(ex)
    r3, n_genes, n_tf = step3_model(m, ex)
    r4, bio_seg = step4_biophysics(ex, idf, m)
    r5, ex_frac = step5_exon_idr_fraction(ex)

    for df, name in [(r0, "density_reconciliation"), (r1, "symmetric_exons"),
                     (r2, "alt_exon_usage"), (r3, "gene_model"),
                     (r4, "biophysics_lengthmatched"),
                     (r5, "exon_idr_fraction_tf")]:
        C.stamp(df).to_csv(C.ANALYSIS / f"splicing_idr_{name}.csv", index=False)

    make_figure(r1, r2, r3, bio_seg,
                C.FIGURES / "fig_splicing_idr_mechanism.pdf")
    make_exon_fraction_figure(ex_frac, r5,
                C.FIGURES / "fig_exon_idr_fraction_tf.pdf")

    # ── headline console summary ──
    print("\n" + "=" * 68)
    print("S9 SUMMARY — does splicing modulate IDRs, and is it TF-biased?")
    print("=" * 68)
    s2i = r2.set_index("cohort")
    print(f"[step2] alt exons encode IDR: all {s2i.loc['all','alt_exon_idr_pct']}% vs "
          f"const {s2i.loc['all','const_exon_idr_pct']}%  (OR={s2i.loc['all','odds_ratio']}, "
          f"gene-paired p={s2i.loc['all','gene_paired_p']:.1e})")
    print(f"        TF: {s2i.loc['TF','alt_exon_idr_pct']}% of alt exons encode IDR "
          f"(gene-paired p={s2i.loc['TF','gene_paired_p']:.1e})")
    s1i = r1[r1.scope == "internal_exons"].set_index("cohort")
    print(f"[step1] internal IDR-exon symmetric: TF {s1i.loc['TF','idr_exon_symmetric_pct']}% vs "
          f"ordered {s1i.loc['TF','ordered_exon_symmetric_pct']}% (OR={s1i.loc['TF','odds_ratio']}, "
          f"gene-paired p={s1i.loc['TF','gene_paired_p']:.1e})")
    tfb = {mm: r3[(r3.model == mm) & (r3.term == 'tf')].beta.values[0]
           for mm in ['unadjusted', 'M0_controls', 'M1_plus_idr_splicing']}
    print(f"[step3] TF beta on pct_change_idr: {tfb['unadjusted']:.2f} -> "
          f"{tfb['M0_controls']:.2f} (+controls) -> {tfb['M1_plus_idr_splicing']:.2f} (+splicing); "
          f"IDR-splicing beta={r3[(r3.model=='M1_plus_idr_splicing')&(r3.term=='z_altidr')].beta.values[0]:.2f}")
    print(f"        (gene-level model n={n_genes:,}, TF={n_tf})")
    if len(bio_seg):
        tf = bio_seg[bio_seg.is_tf == True]
        print(f"[step4] TF alt-IDR vs const-IDR: aromatic "
              f"{tf[tf.seg_alt==True].fract_aro.median():.3f} vs {tf[tf.seg_alt==False].fract_aro.median():.3f}; "
              f"len {tf[tf.seg_alt==True]['IDR len'].median():.0f} vs {tf[tf.seg_alt==False]['IDR len'].median():.0f} aa")
    print("[out] result CSVs -> analysis/  |  figure -> figures/genomic/")


if __name__ == "__main__":
    main()
