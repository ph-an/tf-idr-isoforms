"""
02 (part ①③) — isoform→transcript resolution by CDS structure, + exon region class.

Fills the gap that s3's "best-validating pick" leaves open. When a UniProt isoform
has several sequence-valid ENSTs, we do NOT choose one arbitrarily: we build a
coding-structure signature for each (chrom, strand, ordered CDS-exon intervals,
phases) and decide whether they are the SAME protein-coding structure that differs
only in UTRs (→ collapse) or GENUINELY different CDS exon structures (→ flag
ambiguous, keep all).

resolved_mapping_status ∈
    resolved_unique                 one sequence-valid ENST
    resolved_collapsed_equivalent   >1 valid ENST, identical CDS structure (UTR-only diff)
    ambiguous_multiple_cds_structures  >1 valid ENST, different CDS structures (don't pick)
    unresolved_no_valid_transcript  candidate ENST(s) but none validate
    gene_only_mapping               ENSG known, no isoform-specific ENST
    no_mapping                      no Ensembl mapping

Reuses: s3.load_gencode_translations / s3.match_type / s3.MATCH_RANK,
        isoform_transcript_candidates_long.csv, isoform_transcript_map.csv,
        IDRisoforms_df_geneage.csv, ensembl_transcript_exon_coordinates.parquet,
        tf_idr_isoform_master.parquet, analysis/exon_idr_annotation.csv (s9).

Outputs -> data/genomic/isoform_mapping/  and  figures/genomic/
"""
import sys, hashlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C
from tfidr_pipeline.stages import s3_validate_translations as s3

OUTDIR = C.ISOFORM_MAPPING; OUTDIR.mkdir(exist_ok=True)
FIGDIR = C.FIGURES; FIGDIR.mkdir(parents=True, exist_ok=True)

VALID = C.VALID_SEQUENCE_MATCH_TYPES
OKABE = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "gray": "#7A7A7A"}


def cds_signature(enst, exon_by_enst):
    """(chrom, strand, ((cds_start,cds_end)...)) over coding exons in rank order.
    UTR-only-different transcripts collapse to the same signature; different CDS
    exon structures do not. Returns None if the ENST has no coding exons here."""
    g = exon_by_enst.get(enst)
    if g is None:
        return None
    g = g[g["is_coding_exon"]].sort_values("exon_rank")
    if g.empty:
        return None
    intervals = tuple((int(s), int(e)) for s, e in
                      zip(g["cds_genomic_start"], g["cds_genomic_end"]))
    return (g["chromosome"].iloc[0], g["strand"].iloc[0], intervals)


def sig_hash(sig):
    return None if sig is None else hashlib.sha1(repr(sig).encode()).hexdigest()[:10]


def resolve(filt=None):
    print("[load] inputs ...")
    m = pd.read_csv(C.ISOFORM_MAP)                       # one row / isoform + mapping_status
    cand = pd.read_csv(C.ISOFORM_CANDIDATES)             # one row / candidate ENST
    uni_seq = (pd.read_csv(C.ISOFORM_ANNOTATED)
               .set_index("isoform_accession")["Sequence"].to_dict())
    gc = s3.load_gencode_translations()
    ex = pd.read_parquet(C.INTERIM / "ensembl_transcript_exon_coordinates.parquet")
    exon_by_enst = {e: g for e, g in ex.groupby("ensembl_transcript_id")}
    master = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    meta = master.drop_duplicates("isoform_accession").set_index("isoform_accession")

    cand["ENST_u"] = cand["ENST"].astype(str).str.split(".").str[0]
    cand_by_iso = {iso: list(zip(g["ENST_u"], g["ENSP"], g["ENSG"]))
                   for iso, g in cand.groupby("uniprot_isoform")}

    rows = []
    for r in m.itertuples():
        iso, base, status = r.uniprot_isoform, r.base_accession, r.mapping_status
        gene = meta["gene_name"].get(iso)
        is_tf = bool(meta["is_tf"].get(iso, False))
        tf_family = meta["tf_family"].get(iso) if "tf_family" in meta else None
        ensg = r.ENSG
        rec = dict(base_accession=base, isoform_accession=iso, gene_name=gene,
                   is_tf=is_tf, tf_family=tf_family, ENSG=ensg,
                   ENST_resolved=None, ENSP_resolved=None,
                   n_candidate_ENST=0, n_sequence_valid_ENST=0,
                   n_distinct_cds_structures=0, coding_structure_signature=None,
                   resolved_mapping_status=status, notes="")

        if status in ("gene_only_mapping", "no_mapping"):
            rows.append(rec); continue

        uni = str(uni_seq.get(iso, "") or "")
        cands = cand_by_iso.get(iso, [])
        rec["n_candidate_ENST"] = len({e for e, _, _ in cands})

        # validate every candidate; keep the sequence-valid ones with their score
        valid = []
        for enst, ensp, cg in cands:
            mt, _ = s3.match_type(uni, gc.get(enst, ""))
            if mt in VALID:
                valid.append((enst, ensp, cg, s3.MATCH_RANK[mt]))
        rec["n_sequence_valid_ENST"] = len({e for e, _, _, _ in valid})

        if not valid:
            rec["resolved_mapping_status"] = "unresolved_no_valid_transcript"
            rec["notes"] = "candidate ENST(s) present but none translate to the UniProt isoform"
            rows.append(rec); continue

        # signatures of the valid ENSTs
        sigs = {}
        for enst, ensp, cg, rank in valid:
            sigs.setdefault(enst, (cds_signature(enst, exon_by_enst), ensp, cg, rank))
        n_valid_enst = len(sigs)
        with_sig = {e: v for e, v in sigs.items() if v[0] is not None}
        distinct = {sig_hash(v[0]) for v in with_sig.values()}
        rec["n_distinct_cds_structures"] = len(distinct)
        # representative = best-validating ENST that HAS CDS coordinates, else any
        pool = with_sig or sigs
        rep_enst, (rep_sig, rep_ensp, rep_ensg, _) = sorted(
            pool.items(), key=lambda kv: (-kv[1][3], kv[0]))[0]
        rec["coding_structure_signature"] = sig_hash(rep_sig)
        rec["ENSG"] = rep_ensg or ensg

        if len(distinct) >= 2:
            rec["resolved_mapping_status"] = "ambiguous_multiple_cds_structures"
            rec["notes"] = (f"{n_valid_enst} valid ENSTs, {len(distinct)} distinct CDS "
                            f"structures — NOT auto-picked: {','.join(sorted(sigs))}")
        elif len(distinct) == 1:
            rec["ENST_resolved"], rec["ENSP_resolved"] = rep_enst, rep_ensp
            if n_valid_enst == 1:
                rec["resolved_mapping_status"] = "resolved_unique"
            else:
                rec["resolved_mapping_status"] = "resolved_collapsed_equivalent"
                rec["notes"] = (f"{n_valid_enst} valid ENSTs, identical CDS structure "
                                f"(UTR-only differences): {','.join(sorted(sigs))}")
        else:  # no valid ENST is present in the parsed GENCODE annotation → no coords
            rec["resolved_mapping_status"] = "valid_no_cds_coordinates"
            rec["ENST_resolved"], rec["ENSP_resolved"] = rep_enst, rep_ensp
            rec["notes"] = (f"{n_valid_enst} sequence-valid ENST(s), none in the parsed "
                            f"GENCODE v{C.GENCODE_RELEASE} annotation: {','.join(sorted(sigs))}")
        rows.append(rec)

    res = C.stamp(pd.DataFrame(rows))
    if filt:
        keep = set(filt)
        res = res[res.isoform_accession.isin(keep) | res.base_accession.isin(keep)
                  | res.gene_name.isin(keep)]
    res.to_parquet(OUTDIR / "isoform_transcript_mapping_resolved.parquet", index=False)
    res.to_csv(OUTDIR / "isoform_transcript_mapping_resolved.csv", index=False)
    return res


def summaries(res):
    order = ["resolved_unique", "resolved_collapsed_equivalent",
             "ambiguous_multiple_cds_structures", "valid_no_cds_coordinates",
             "unresolved_no_valid_transcript", "gene_only_mapping", "no_mapping"]
    by = (res.resolved_mapping_status.value_counts()
          .reindex(order).fillna(0).astype(int).rename("n").reset_index())
    by.columns = ["resolved_mapping_status", "n"]
    by["pct"] = (by.n / len(res) * 100).round(1)
    by.to_csv(OUTDIR / "resolved_mapping_status_summary.csv", index=False)

    ct = (pd.crosstab(res.resolved_mapping_status, res.is_tf)
          .reindex(order).fillna(0).astype(int))
    ct.columns = ["non_TF" if c is False else "TF" for c in ct.columns]
    ct.to_csv(OUTDIR / "resolved_mapping_status_by_tf.csv")

    MAIN = {"resolved_unique", "resolved_collapsed_equivalent"}
    res["in_main_analysis"] = res.resolved_mapping_status.isin(MAIN)
    ret = res.groupby("is_tf")["in_main_analysis"].agg(["sum", "count"])
    ret["pct_retained"] = (ret["sum"] / ret["count"] * 100).round(1)
    return by, ct, ret


def make_qc_figure(ct):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[fig] matplotlib unavailable ({e})"); return
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    order = list(ct.index)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    y = np.arange(len(order)); h = 0.38
    ax.barh(y + h/2, ct.get("TF", pd.Series(0, index=order)).values, h,
            color=OKABE["verm"], label="TF")
    ax.barh(y - h/2, ct.get("non_TF", pd.Series(0, index=order)).values, h,
            color=OKABE["blue"], label="non-TF")
    for i, s in enumerate(order):
        tf = int(ct.get("TF", pd.Series(0, index=order)).get(s, 0))
        nt = int(ct.get("non_TF", pd.Series(0, index=order)).get(s, 0))
        ax.text(tf + 60, i + h/2, f"{tf:,}", va="center", fontsize=7)
        ax.text(nt + 60, i - h/2, f"{nt:,}", va="center", fontsize=7)
    ax.set_yticks(y); ax.set_yticklabels([s.replace("_", " ") for s in order])
    ax.set_xlabel("number of UniProt isoforms"); ax.invert_yaxis()
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Isoform→transcript resolution by CDS structure", loc="left")
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig_mapping_resolution_qc.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] wrote fig_mapping_resolution_qc.pdf")


def exon_region_class():
    """Part ③ — finish the exon→protein map with a region class + genomic coords."""
    ann_path = C.ANALYSIS / "exon_idr_annotation.csv"
    if not ann_path.exists():
        print("[③] exon_idr_annotation.csv not found (run s9 first); skipping"); return None
    ex = pd.read_csv(ann_path)
    # region_class from the IDR fraction of the exon's residues
    ex["region_class"] = np.select(
        [ex.idr_frac >= 0.99, ex.idr_frac <= 0.01],
        ["entirely_IDR", "entirely_nonIDR"], default="crosses_IDR_boundary")
    # attach genomic coordinates from the CDS→protein map
    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")
    cds = cds[["ensembl_transcript_id", "exon_id", "chromosome", "strand",
               "cds_genomic_start", "cds_genomic_end"]].drop_duplicates(["ensembl_transcript_id", "exon_id"])
    out = ex.merge(cds, left_on=["ENST", "exon_id"],
                   right_on=["ensembl_transcript_id", "exon_id"], how="left")
    out.to_parquet(OUTDIR / "exon_region_class.parquet", index=False)
    vc = out.region_class.value_counts()
    print(f"[③] exon_region_class.parquet  {len(out):,} exons  "
          f"({dict(vc)})")
    return out


def main(filt=None):
    res = resolve(filt)
    by, ct, ret = summaries(res)
    make_qc_figure(ct)
    exon_region_class()

    print("\n" + "=" * 66)
    print("① ISOFORM→TRANSCRIPT RESOLUTION (by CDS structure)")
    print("=" * 66)
    print(by.to_string(index=False))
    n_main = int(res.resolved_mapping_status.isin(
        ["resolved_unique", "resolved_collapsed_equivalent"]).sum())
    print(f"\nMAIN analysis set (unique + collapsed-equivalent): {n_main:,} "
          f"({n_main/len(res)*100:.1f}%)")
    amb = int((res.resolved_mapping_status == "ambiguous_multiple_cds_structures").sum())
    coll = int((res.resolved_mapping_status == "resolved_collapsed_equivalent").sum())
    print(f"collapsed UTR-only-equivalent (rescued from ambiguity): {coll:,}")
    print(f"genuinely ambiguous (multiple CDS structures, NOT picked): {amb:,}")
    print("\nretention into MAIN by TF status:")
    print(ret.to_string())
    print(f"\n[out] {OUTDIR}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
