"""
Stage D (part 1) — fetch the GENCODE annotation once.

Downloads three modest files (no genome FASTA, no 16k REST calls needed):
  * <release>.annotation.gtf.gz         -> exon/CDS genomic coordinates (comprehensive)
  * <release>.pc_translations.fa.gz     -> per-ENST protein sequences (Stage C
                                            at-scale validation vs UniProt)
  * <release>.pc_transcripts.fa.gz      -> per-ENST spliced transcript nucleotides
                                            + CDS offsets (s4b DNA round-trip that
                                            self-verifies the genomic projection)

Idempotent: skips anything already cached.
"""
import sys, urllib.request, time
sys.stdout.reconfigure(encoding="utf-8")
from tfidr_pipeline import config as C


def fetch(url, dest):
    if dest.exists():
        print(f"[cache] {dest.name} ({dest.stat().st_size/1e6:.0f} MB)")
        return
    print(f"[download] {url}")
    t0 = time.time()
    urllib.request.urlretrieve(url, dest)
    print(f"[download] {dest.name} -> {dest.stat().st_size/1e6:.0f} MB in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    fetch(C.GENCODE_GTF_URL, C.GENCODE_GTF_GZ)
    fetch(C.GENCODE_TRANSLATIONS_URL, C.GENCODE_TRANSLATIONS_GZ)
    fetch(C.GENCODE_PC_TRANSCRIPTS_URL, C.GENCODE_PC_TRANSCRIPTS_GZ)
    print("Done.")
