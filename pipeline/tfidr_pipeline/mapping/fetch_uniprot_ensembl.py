"""
STEP 1-2 of the UniProt-isoform -> Ensembl-transcript workflow.

Bulk-downloads the isoform-resolved Ensembl cross-references for the entire
human SwissProt proteome from UniProt, then parses them into a tidy table:

    uniprot_isoform   ENST   ENSP   ENSG

UniProt's `xref_ensembl_full` TSV field encodes each cross-reference as
    "ENST00000269305.9; ENSP00000269305.4; ENSG00000141510.20. [P04637-1]"
The bracketed [P04637-1] is the *isoform-specific* tag — this is what lets us
attach a genomic transcript to a particular UniProt isoform rather than just
to the gene.

Output: data/raw/uniprot/uniprot_ensembl_xref_raw.tsv   (raw download, so we never refetch)
        cache/uniprot_ensembl_isoform_pairs.csv  (parsed long table)
"""
import sys, io, re, time, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import requests
import pandas as pd
from tfidr_pipeline import config as C

CACHE  = C.CACHE
RAW_TSV   = C.UNIPROT_XREF_TSV
PAIRS_CSV = CACHE / "uniprot_ensembl_isoform_pairs.csv"

STREAM = "https://rest.uniprot.org/uniprotkb/stream"


def fetch_raw():
    """Stream the whole reviewed human proteome's Ensembl xref field once."""
    if RAW_TSV.exists():
        print(f"[cache] using existing {RAW_TSV.name} ({RAW_TSV.stat().st_size/1e6:.1f} MB)")
        return RAW_TSV.read_text(encoding="utf-8")

    params = {
        "query":  "reviewed:true AND organism_id:9606",
        "fields": "accession,xref_ensembl_full",
        "format": "tsv",
    }
    print("[fetch] streaming Ensembl cross-references for reviewed human proteome ...")
    t0 = time.time()
    with requests.get(STREAM, params=params, timeout=300, stream=True) as r:
        r.raise_for_status()
        text = r.text
    RAW_TSV.write_text(text, encoding="utf-8")
    n = text.count("\n")
    print(f"[fetch] done in {time.time()-t0:.1f}s  |  {n:,} entry lines  |  saved {RAW_TSV.name}")
    return text


# One xref token: "ENST...; ENSP...; ENSG.... [ISOFORM]"
_TOKEN = re.compile(
    r"(ENST\d+(?:\.\d+)?)\s*;\s*"
    r"(ENSP\d+(?:\.\d+)?)\s*;\s*"
    r"(ENSG\d+(?:\.\d+)?)\.?"
    r"(?:\s*\[([A-Z0-9]+-\d+)\])?"
)


def parse(text):
    """Turn the raw TSV into a long table of isoform<->transcript pairs."""
    df = pd.read_csv(io.StringIO(text), sep="\t")
    df.columns = ["Entry", "Ensembl"]
    df = df[df["Ensembl"].notna()]

    rows = []
    for entry, cell in zip(df["Entry"], df["Ensembl"]):
        for m in _TOKEN.finditer(cell):
            enst, ensp, ensg, iso = m.groups()
            rows.append({
                "base_accession":  entry,
                "isoform_tag":     iso if iso else entry,  # unbracketed => single isoform
                "ENST":            enst,
                "ENSP":            ensp,
                "ENSG":            ensg,
            })
    pairs = pd.DataFrame(rows)
    # strip version suffixes into separate unversioned columns (join-friendly)
    for col in ["ENST", "ENSP", "ENSG"]:
        pairs[col + "_nover"] = pairs[col].str.split(".").str[0]
    pairs.to_csv(PAIRS_CSV, index=False)
    return pairs


if __name__ == "__main__":
    text  = fetch_raw()
    pairs = parse(text)
    print(f"\nParsed {len(pairs):,} isoform<->transcript pairs")
    print(f"  unique base accessions : {pairs['base_accession'].nunique():,}")
    print(f"  unique isoform tags    : {pairs['isoform_tag'].nunique():,}")
    print(f"  unique ENST            : {pairs['ENST'].nunique():,}")
    print(f"\nSaved -> {PAIRS_CSV}")
    print(pairs.head(12).to_string(index=False))
