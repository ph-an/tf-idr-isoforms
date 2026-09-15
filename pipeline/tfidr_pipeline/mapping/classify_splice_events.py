"""
02 (part ②) — canonical↔alternative splice-event classification + IDR consequence.

For every gene whose UniProt canonical isoform AND ≥1 alternative isoform both have
a MAIN-resolved transcript (from resolve_isoform_transcripts.py), compare their
coding-exon structures, classify the splice event(s) that distinguish them, and ask
what each event does to the isoform's IDRs.

Event taxonomy (v1): exon_skip · alt_5p_splice_site · alt_3p_splice_site ·
alt_first_exon · alt_last_exon · intron_retention · complex.

IDR consequence per differential exon (idr_frac = fraction of its residues that are
IDR): removes_IDR / adds_IDR / shortens_IDR / extends_IDR / alters_ordered.

The biological readout is P(IDR-altering | event type), split by TF / non-TF.

Reuses: isoform_mapping/isoform_transcript_mapping_resolved.parquet,
        transcript_cds_to_protein_coordinate_map.parquet, MetaPredict IDRome,
        tf_idr_isoform_master.parquet.

Outputs -> data/genomic/isoform_mapping/  and  figures/genomic/
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from tfidr_pipeline import config as C
from tfidr_pipeline import utils as U

OUTDIR = C.ISOFORM_MAPPING; OUTDIR.mkdir(exist_ok=True)
FIGDIR = C.FIGURES; FIGDIR.mkdir(parents=True, exist_ok=True)
MAIN = {"resolved_unique", "resolved_collapsed_equivalent"}
OKABE = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "gray": "#7A7A7A"}


def five_p(exon, strand):   # 5' (upstream) boundary in transcription direction
    return exon[0] if strand == "+" else exon[1]


def three_p(exon, strand):  # 3' (downstream) boundary
    return exon[1] if strand == "+" else exon[0]


def idr_overlap_frac(aa_start, aa_end, segs):
    L = aa_end - aa_start + 1
    return U.overlap_length(aa_start, aa_end, segs) / L if L > 0 and segs else 0.0


def classify_pair(C_ex, A_ex, strand):
    """C_ex, A_ex: lists of dicts {iv:(start,end), aa:(s,e)} sorted 5'->3'.
    Returns list of events: (event_type, exon_iv, exon_aa, present_in)."""
    Civ = [e["iv"] for e in C_ex]; Aiv = [e["iv"] for e in A_ex]
    Cset, Aset = set(Civ), set(Aiv)
    c_first, c_last = (Civ[0], Civ[-1]) if Civ else (None, None)
    a_first, a_last = (Aiv[0], Aiv[-1]) if Aiv else (None, None)
    c_only = [e for e in C_ex if e["iv"] not in Aset]
    a_only = [e for e in A_ex if e["iv"] not in Cset]
    events = []

    def contains(iv, others):
        return [o for o in others if o[0] >= iv[0] and o[1] <= iv[1] and o != iv]

    # 0) intron retention — one exon fully spans ≥2 exons of the other transcript
    consumed_c, consumed_a = set(), set()
    for e in a_only:                              # alternative retains an intron
        sp = contains(e["iv"], Civ)
        if len(sp) >= 2:
            events.append(("intron_retention", e["iv"], e["aa"], "alternative"))
            consumed_a.add(e["iv"]); consumed_c.update(sp)
    for e in c_only:                              # canonical is the intron-retaining form
        sp = contains(e["iv"], Aiv)
        if len(sp) >= 2:
            events.append(("intron_retention", e["iv"], e["aa"], "canonical"))
            consumed_c.add(e["iv"]); consumed_a.update(sp)
    c_only = [e for e in c_only if e["iv"] not in consumed_c]
    a_only = [e for e in a_only if e["iv"] not in consumed_a]

    # 1) overlapping differential exons that share one boundary = alt splice site
    used_a, remaining_c = set(), []
    for e in c_only:
        c = e["iv"]
        match = next((a for a in a_only if a["iv"] not in used_a and U.intervals_overlap(c, a["iv"])), None)
        if match is None:
            remaining_c.append(e); continue
        a = match["iv"]; used_a.add(a)
        if five_p(c, strand) == five_p(a, strand) and three_p(c, strand) != three_p(a, strand):
            et = "alt_5p_splice_site"            # donor varies
        elif three_p(c, strand) == three_p(a, strand) and five_p(c, strand) != five_p(a, strand):
            et = "alt_3p_splice_site"            # acceptor varies
        else:
            et = "complex"
        events.append((et, c, e["aa"], "both"))

    # 2) non-overlapping leftovers → skip / terminal
    def leftover(exs, first, last, present):
        for e in exs:
            iv = e["iv"]
            et = ("alt_first_exon" if iv == first else
                  "alt_last_exon" if iv == last else "exon_skip")
            events.append((et, iv, e["aa"], present))

    leftover(remaining_c, c_first, c_last, "canonical")
    leftover([a for a in a_only if a["iv"] not in used_a], a_first, a_last, "alternative")
    return events


def main():
    print("[load] inputs ...")
    res = pd.read_parquet(OUTDIR / "isoform_transcript_mapping_resolved.parquet")
    res = res[res.resolved_mapping_status.isin(MAIN) & res.ENST_resolved.notna()].copy()
    res["ENST"] = res.ENST_resolved.astype(str).str.split(".").str[0]

    master = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    canon = master.drop_duplicates("isoform_accession").set_index("isoform_accession")
    res["is_canonical"] = res.isoform_accession.map(canon["is_canonical"]).fillna(False)

    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")
    cds = cds.sort_values(["ensembl_transcript_id", "exon_rank"])
    strand_of, cds_exons = {}, {}
    for enst, g in cds.groupby("ensembl_transcript_id"):
        strand_of[enst] = g.strand.iloc[0]
        cds_exons[enst] = [{"iv": (int(s), int(e)), "aa": (int(a0), int(a1))}
                           for s, e, a0, a1 in zip(g.cds_genomic_start, g.cds_genomic_end,
                                                   g.aa_start, g.aa_end)]

    idr_by = U.idr_segments_by_isoform()

    rows = []
    n_genes = n_pairs = 0
    for gene, g in res.groupby("base_accession"):
        canon_rows = g[g.is_canonical]
        alts = g[~g.is_canonical]
        if canon_rows.empty or alts.empty:
            continue
        cr = canon_rows.iloc[0]
        c_iso, c_enst = cr.isoform_accession, cr.ENST
        if c_enst not in cds_exons:
            continue
        n_genes += 1
        C_ex = cds_exons[c_enst]; strand = strand_of.get(c_enst, "+")
        for a in alts.itertuples():
            if a.ENST not in cds_exons or a.ENST == c_enst:
                continue
            n_pairs += 1
            for et, iv, aa, present in classify_pair(C_ex, cds_exons[a.ENST], strand):
                owner_iso = c_iso if present in ("canonical", "both") else a.isoform_accession
                frac = idr_overlap_frac(aa[0], aa[1], idr_by.get(owner_iso, []))
                is_idr_exon = frac >= 0.5
                if present == "canonical":
                    cons = "removes_IDR" if is_idr_exon else "removes_ordered"
                elif present == "alternative":
                    cons = "adds_IDR" if is_idr_exon else "adds_ordered"
                else:  # alt splice site: shorten vs extend by exon length
                    cons = ("resizes_IDR" if is_idr_exon else "resizes_ordered")
                rows.append(dict(
                    base_accession=gene, gene_name=cr.gene_name, is_tf=bool(cr.is_tf),
                    canonical_isoform=c_iso, alternative_isoform=a.isoform_accession,
                    canonical_ENST=c_enst, alternative_ENST=a.ENST,
                    event_type=et, exon_genomic_start=iv[0], exon_genomic_end=iv[1],
                    exon_aa_start=aa[0], exon_aa_end=aa[1], present_in=present,
                    exon_idr_frac=round(frac, 3), is_idr_exon=is_idr_exon,
                    idr_consequence=cons))

    ev = C.stamp(pd.DataFrame(rows))
    ev.to_parquet(OUTDIR / "splice_events.parquet", index=False)
    ev.to_csv(OUTDIR / "splice_events.csv", index=False)

    # ── P(IDR-altering | event type), TF vs non-TF ──
    summ = []
    for coh, sub in [("all", ev), ("TF", ev[ev.is_tf]), ("non_TF", ev[~ev.is_tf])]:
        for et, g in sub.groupby("event_type"):
            summ.append(dict(cohort=coh, event_type=et, n_events=len(g),
                             pct_idr_altering=round(g.is_idr_exon.mean() * 100, 1)))
    summ = pd.DataFrame(summ)
    C.stamp(summ).to_csv(OUTDIR / "splice_event_idr_summary.csv", index=False)

    # TF-enrichment of IDR-altering per event type (fisher, TF vs non-TF)
    make_figure(ev, summ)

    print(f"\n[out] genes compared: {n_genes:,} | canonical↔alt pairs: {n_pairs:,} "
          f"| events: {len(ev):,}")
    print("\nP(IDR-altering | event type):")
    piv = summ.pivot(index="event_type", columns="cohort", values="pct_idr_altering")
    cnt = ev.event_type.value_counts()
    piv["n_all"] = cnt
    print(piv.sort_values("n_all", ascending=False).to_string())
    print("\nIDR consequence breakdown (all):")
    print(ev.idr_consequence.value_counts().to_string())
    print(f"\n[out] {OUTDIR}")


def make_figure(ev, summ):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[fig] matplotlib unavailable ({e})"); return
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    order = (ev.event_type.value_counts().index.tolist())
    piv = summ.pivot(index="event_type", columns="cohort", values="pct_idr_altering").reindex(order)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    y = np.arange(len(order)); h = 0.38
    ax.barh(y + h/2, piv["TF"].values, h, color=OKABE["verm"], label="TF")
    ax.barh(y - h/2, piv["non_TF"].values, h, color=OKABE["blue"], label="non-TF")
    for i, et in enumerate(order):
        n = int((ev.event_type == et).sum())
        ax.text(1, i + 0.28, f"n={n:,}", fontsize=6.5, va="center", color="#444")
    ax.axvline(50, ls=":", lw=1, color="k", alpha=.4)
    ax.set_yticks(y); ax.set_yticklabels([e.replace("_", " ") for e in order])
    ax.set_xlabel("% of events that are IDR-altering"); ax.set_xlim(0, 100)
    ax.invert_yaxis(); ax.legend(frameon=False, loc="lower right")
    ax.set_title("Splice-event type → IDR consequence (canonical vs alternative)", loc="left")
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig_splice_event_idr.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[fig] wrote fig_splice_event_idr.pdf")


if __name__ == "__main__":
    main()
