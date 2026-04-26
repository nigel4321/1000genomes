"""Truth-set BED tracks (M11).

Two BED4 tracks are emitted alongside each person VCF:

* `<truth_dir>/person_NNNN.golden.bed` — the curated set of variants
  the spec calls "golden truth". Includes the per-person highlighted
  ClinVar variant, every cohort row carrying a ClinVar / dbSNP /
  COSMIC annotation, and every structural variant. This is the ground
  truth a downstream caller would aim to recover.
* `<truth_dir>/person_NNNN.noise.bed` — every per-call perturbation
  introduced by the M9 sequencing-error model: GT flips and
  dropouts. Each line records the truth GT and the called GT so a
  caller's accuracy can be graded against the model's known noise.

BED4 is the lowest-common-denominator BED format that downstream
tools accept (`bedtools`, IGV, UCSC). The 4th column is a
semicolon-separated `key=value` payload so the tag set can grow
without breaking parsers.
"""

from __future__ import annotations

from pathlib import Path


# Golden record categories — every variant lands at most one tag, in
# this priority order (highlighted is rarest, SV the most common).
GOLDEN_CATEGORIES = ("HIGHLIGHTED", "CLINVAR", "COSMIC", "RSID", "SV")


def classify_golden(variant: dict, is_hi: bool) -> str | None:
    """Return the golden category for a record, or None if it's not a
    golden-set member.

    Tag priority: HIGHLIGHTED > CLINVAR > COSMIC > RSID > SV.
    A row that has both CLINVAR and RSID (e.g. an injected ClinVar
    record whose ID is `rs…`) gets tagged CLINVAR — the more specific
    annotation wins.
    """
    if is_hi:
        return "HIGHLIGHTED"
    if variant.get("clnsig") and variant["clnsig"] != ".":
        return "CLINVAR"
    if variant.get("cosmic_id"):
        return "COSMIC"
    if variant.get("svtype"):
        return "SV"
    rid = variant.get("id") or ""
    if rid.startswith("rs"):
        return "RSID"
    return None


def _bed_end(variant: dict) -> int:
    """Half-open end coordinate for the variant.

    SV records expose `end` already (POS + |SVLEN|); SNV/indel records
    use `(pos - 1) + len(ref)` so the BED span equals len(ref) under
    the standard 0-based half-open convention (a 1-base SNV at pos 1000
    becomes [999, 1000)).
    """
    if variant.get("end") is not None:
        return int(variant["end"])
    return (int(variant["pos"]) - 1) + max(1, len(variant.get("ref", "N")))


def _format_payload(items) -> str:
    """`[(key, value), ...]` → `key1=value1;key2=value2`. Skips `None`
    values and replaces tab / semicolon / newline so the BED stays
    parseable."""
    out = []
    for k, v in items:
        if v is None:
            continue
        s = str(v).replace("\t", " ").replace(";", ",") \
            .replace("\n", " ")
        out.append(f"{k}={s}")
    return ";".join(out) if out else "."


def golden_bed_line(variant: dict, category: str, gt: str) -> str:
    """One BED4 row tagging a golden-set record."""
    chrom = variant["chrom"]
    start = int(variant["pos"]) - 1  # BED is 0-based half-open
    end = _bed_end(variant)
    payload = _format_payload([
        ("flag", category),
        ("id", variant.get("id") or "."),
        ("ref", variant.get("ref")),
        ("alt", ",".join(variant.get("alts", []))),
        ("gt", gt),
        ("clnsig", variant.get("clnsig")),
        ("clndn", variant.get("clndn")),
        ("cosmic_id", variant.get("cosmic_id")),
        ("cosmic_gene", variant.get("cosmic_gene")),
        ("svtype", variant.get("svtype")),
        ("svlen", variant.get("svlen")),
    ])
    return f"{chrom}\t{start}\t{end}\t{payload}"


def noise_bed_line(variant: dict, kind: str, truth_gt: str,
                   called_gt: str) -> str:
    """One BED4 row tagging a noise event (flip or dropout).

    `kind` ∈ {"FLIP", "DROPOUT"}.
    """
    chrom = variant["chrom"]
    start = int(variant["pos"]) - 1
    end = _bed_end(variant)
    payload = _format_payload([
        ("flag", kind),
        ("id", variant.get("id") or "."),
        ("ref", variant.get("ref")),
        ("alt", ",".join(variant.get("alts", []))),
        ("truth_gt", truth_gt),
        ("called_gt", called_gt),
    ])
    return f"{chrom}\t{start}\t{end}\t{payload}"


class TruthBedWriter:
    """Manages the pair of golden / noise BED files for one person.

    Lines are buffered in memory and flushed in genomic order on
    `close()` so the BED is `sort -k1,1 -k2,2n`-friendly without a
    follow-up shell sort.
    """

    def __init__(self, golden_path: Path, noise_path: Path,
                 contig_order: dict | None = None):
        self.golden_path = Path(golden_path)
        self.noise_path = Path(noise_path)
        self.contig_order = contig_order or {}
        self._golden: list[tuple[tuple, str]] = []
        self._noise: list[tuple[tuple, str]] = []
        self.golden_count = 0
        self.noise_count = 0

    def _key(self, chrom: str, start: int) -> tuple:
        return (self.contig_order.get(chrom, len(self.contig_order)),
                chrom, start)

    def add_golden(self, variant: dict, category: str,
                   gt: str) -> None:
        line = golden_bed_line(variant, category, gt)
        start = int(variant["pos"]) - 1
        self._golden.append((self._key(variant["chrom"], start), line))
        self.golden_count += 1

    def add_noise(self, variant: dict, kind: str, truth_gt: str,
                  called_gt: str) -> None:
        line = noise_bed_line(variant, kind, truth_gt, called_gt)
        start = int(variant["pos"]) - 1
        self._noise.append((self._key(variant["chrom"], start), line))
        self.noise_count += 1

    def _flush(self, path: Path, rows: list) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        rows.sort(key=lambda r: r[0])
        with open(path, "w") as fh:
            for _, line in rows:
                fh.write(line + "\n")

    def close(self) -> None:
        self._flush(self.golden_path, self._golden)
        self._flush(self.noise_path, self._noise)

    # Allow `with TruthBedWriter(...) as w:` ergonomics.
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
