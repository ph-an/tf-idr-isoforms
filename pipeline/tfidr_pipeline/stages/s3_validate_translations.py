"""
Stage C at scale — validate each mapped ENST's protein against the UniProt
isoform sequence, using the local GENCODE pc_translations FASTA (no REST calls).

For every isoform in isoform_transcript_map.csv that carries an ENST, we compare
the GENCODE translation to the UniProt isoform sequence we already hold in
IDRisoforms_df_geneage.csv and assign a sequence_match_type:

    exact              identical sequences
    near_exact         >=99% identity (equal length), or one a clean prefix of
                       the other (handles trailing-stop / annotation trims)
    length_match_only  same length but <99% identity
    mismatch           different sequence, no clean relationship
    no_translation     ENST absent from the GENCODE translations FASTA

Only ``exact`` matches are considered valid. Other match classes remain in the
output for QC. We also finish the transcript mapping (Susie's "validation picks the transcript"
principle): for isoforms with several candidate ENSTs we validate ALL of them and
pick the best-validating one as `resolved_transcript`. The original primary-ENST
columns (ENST, sequence_match_type, is_validated_mapping) are preserved unchanged
so existing outputs stay reproducible; downstream genomic projection (s4b) uses
the resolved transcript.

Output: interim/mapping_validation_scaled.parquet
"""
import sys, gzip
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
from tfidr_pipeline import config as C

OUT = C.INTERIM / "mapping_validation_scaled.parquet"

# rank match types best -> worst so we can pick the best-validating candidate
MATCH_RANK = {"exact": 3, "near_exact": 2, "length_match_only": 1,
              "mismatch": 0, "no_translation": -1}


def load_gencode_translations():
    """{ENST_unversioned: protein_seq} from pc_translations FASTA.
    Header: >ENSP..|ENST..|ENSG..|OTT..|OTT..|txname|gene|length"""
    seqs, enst, buf = {}, None, []
    with gzip.open(C.GENCODE_TRANSLATIONS_GZ, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if enst is not None:
                    seqs[enst] = "".join(buf)
                parts = line[1:].split("|")
                enst = next((p.split(".")[0] for p in parts if p.startswith("ENST")), None)
                buf = []
            else:
                buf.append(line.strip())
        if enst is not None:
            seqs[enst] = "".join(buf)
    return seqs


def match_type(uni, ens):
    if not ens:
        return "no_translation", 0.0
    if uni == ens:
        return "exact", 1.0
    # clean prefix relationship (trailing stop or annotation trim)
    if ens.startswith(uni) or uni.startswith(ens):
        return "near_exact", 1.0
    if len(uni) == len(ens):
        ident = sum(a == b for a, b in zip(uni, ens)) / len(uni)
        if ident >= 0.99:
            return "near_exact", ident
        return "length_match_only", ident
    return "mismatch", 0.0


if __name__ == "__main__":
    print("[load] GENCODE translations ...")
    gc = load_gencode_translations()
    print(f"[load] {len(gc):,} ENST translations")

    m = pd.read_csv(C.ISOFORM_MAP)
    uni_seq = (pd.read_csv(C.ISOFORM_ANNOTATED)
               .set_index("isoform_accession")["Sequence"].to_dict())

    # candidate ENSTs per isoform (one row per candidate) -> best-validating pick
    cand = pd.read_csv(C.ISOFORM_CANDIDATES)
    cand_by_iso = {iso: g["ENST"].dropna().tolist()
                   for iso, g in cand.groupby("uniprot_isoform")}

    def best_candidate(iso, uni):
        """Validate every candidate ENST and return the best-scoring one."""
        best = (None, "no_translation", 0.0, -2)  # enst, mtype, ident, rank
        for censt in cand_by_iso.get(iso, []):
            e = str(censt).split(".")[0]
            mt, idt = match_type(uni, gc.get(e, ""))
            rank = MATCH_RANK[mt]
            if (rank, idt) > (best[3], best[2]):
                best = (e, mt, idt, rank)
        return best

    rows = []
    for r in m.itertuples():
        enst = None if pd.isna(r.ENST) else str(r.ENST).split(".")[0]
        uni  = str(uni_seq.get(r.uniprot_isoform, "") or "")
        ens  = gc.get(enst, "") if enst else ""
        mtype, ident = match_type(uni, ens)             # primary-ENST (unchanged)

        # resolve across ALL candidates (Susie: validation selects the transcript)
        r_enst, r_mtype, r_ident, _ = best_candidate(r.uniprot_isoform, uni) if uni else (None, "no_translation", 0.0, -2)
        r_validated = r_mtype in C.VALID_SEQUENCE_MATCH_TYPES

        rows.append({
            "uniprot_isoform": r.uniprot_isoform,
            "base_accession":  r.base_accession,
            "ENST":            r.ENST,
            "mapping_status":  r.mapping_status,
            "uniprot_isoform_length_aa": len(uni) if uni else None,
            "ensembl_translation_length_aa": len(ens) if ens else None,
            "length_difference_aa": (len(ens) - len(uni)) if (uni and ens) else None,
            "sequence_identity": round(ident, 4),
            "sequence_match_type": mtype,
            "is_validated_mapping": mtype in C.VALID_SEQUENCE_MATCH_TYPES,
            # ── resolved transcript (best-validating candidate) ──
            "resolved_transcript":   r_enst,
            "resolved_match_type":   r_mtype,
            "resolved_identity":     round(r_ident, 4),
            "resolved_is_validated": r_validated,
            "resolved_by_validation": bool(r_validated and r_enst is not None
                                           and r_enst != enst),
        })
    val = pd.DataFrame(rows)
    val.to_parquet(OUT, index=False)

    print(f"\n[out] {len(val):,} isoforms -> {OUT.name}\n")
    print("sequence_match_type breakdown:")
    vc = val["sequence_match_type"].value_counts()
    for k, v in vc.items():
        print(f"  {k:<18} {v:>7,}  ({v/len(val)*100:5.1f}%)")

    print("\nvalidated (exact only) among isoforms that HAVE an ENST:")
    has = val[val["ENST"].notna()]
    ok = has["is_validated_mapping"].sum()
    print(f"  {ok:,}/{len(has):,}  ({ok/len(has)*100:.1f}%)")

    print("\nby mapping_status (validated fraction):")
    g = (has.groupby("mapping_status")["is_validated_mapping"]
         .agg(["sum", "count"]))
    g["pct"] = (g["sum"] / g["count"] * 100).round(1)
    print(g.to_string())

    # ── transcript resolution (best-validating candidate across all ENSTs) ──
    prim_ok = int(val["is_validated_mapping"].sum())
    res_ok  = int(val["resolved_is_validated"].sum())
    switched = int(val["resolved_by_validation"].sum())
    print(f"\nresolved-transcript validation: {res_ok:,} isoforms "
          f"(+{res_ok - prim_ok:,} vs primary-ENST {prim_ok:,})")
    print(f"isoforms whose resolved transcript differs from the primary pick: {switched:,}")
