#!/usr/bin/env python3
"""QC + validation for a single VCF before downstream processing.

Two tiers of check:

Hard (errors) — these mean the file is unreadable or unusable downstream:
  - file exists and is non-empty
  - tabix index is present alongside the VCF
  - bcftools can parse the header
  - header has a #CHROM line with sample columns
  - index reports at least one contig
  - file contains at least one variant record

Soft (warnings) — the file is readable but its attributes do not match
what a human-genome VCF typically looks like. These never abort the
pipeline; they surface in qc_report.md for the user to eyeball:
  - reference field names a recognised human build
    (hg19/GRCh37/b37/hs37d5 or hg38/GRCh38/b38)
  - indexed contigs are among the expected human chromosomes
  - FORMAT=GT is declared
  - INFO declares AF, or both AC and AN (required by scan_variant)
  - variant positions fit inside plausible human chromosome lengths

Always writes the JSON output file. In --strict mode, also exits 1 when
any hard check fails, which aborts the Nextflow workflow.
"""

import argparse
import json
import os
import re
import subprocess
import sys


# Canonical (non-prefixed) human chromosome names, plus chr-prefixed forms.
_HUMAN_CANONICAL = set(str(i) for i in range(1, 23)) | {"X", "Y", "M", "MT"}
HUMAN_CHROMS = _HUMAN_CANONICAL | {f"chr{c}" for c in _HUMAN_CANONICAL}

# Lowercased substrings that indicate a recognised human reference build.
HUMAN_REF_TOKENS = {
    "GRCh37": ("grch37", "b37", "hs37d5", "hg19"),
    "GRCh38": ("grch38", "b38", "hg38", "hs38"),
}

# Upper bound for plausible positions on any human chromosome (GRCh37 chr1
# is ~249 Mb — anything past ~300 Mb indicates wrong genome / scaffold).
MAX_HUMAN_CHROM_LENGTH = 300_000_000


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _detect_build(reference: str) -> str | None:
    ref = reference.lower()
    for build, tokens in HUMAN_REF_TOKENS.items():
        if any(tok in ref for tok in tokens):
            return build
    return None


def _parse_header(header_text: str) -> dict:
    chrom_line = None
    info_ids: set[str] = set()
    reference = ""
    has_gt = False
    for ln in header_text.splitlines():
        if ln.startswith("#CHROM"):
            chrom_line = ln
        elif ln.startswith("##reference"):
            reference = ln.split("=", 1)[1].strip() if "=" in ln else ""
        elif ln.startswith("##INFO=<ID="):
            m = re.match(r"##INFO=<ID=([A-Za-z0-9_]+)", ln)
            if m:
                info_ids.add(m.group(1))
        elif ln.startswith("##FORMAT=<ID=GT"):
            has_gt = True
    samples: list[str] = []
    if chrom_line:
        fields = chrom_line.split("\t")
        samples = fields[9:] if len(fields) > 9 else []
    return {
        "chrom_line_present": chrom_line is not None,
        "samples": samples,
        "reference": reference,
        "info_ids": sorted(info_ids),
        "has_format_gt": has_gt,
    }


def _parse_index(idx_stdout: str) -> tuple[list[str], int | None]:
    """Extract (contigs, total_records) from `bcftools index -s` output."""
    contigs = []
    total = 0
    have_counts = False
    for line in idx_stdout.strip().splitlines():
        parts = line.split()
        if not parts:
            continue
        contigs.append(parts[0])
        if len(parts) >= 3:
            try:
                total += int(parts[2])
                have_counts = True
            except ValueError:
                pass
    return contigs, (total if have_counts else None)


def _max_position_per_contig(vcf: str, contigs: list[str]) -> dict[str, int]:
    """Return the largest POS per contig via tabix seek + tail."""
    out: dict[str, int] = {}
    for c in contigs:
        # Stream from the given contig, grab only POS, keep the max. bcftools
        # streams in order so the last non-empty line holds the max.
        r = _run(["bcftools", "query", "-f", "%POS\n", "-r", c, vcf])
        last = ""
        for line in r.stdout.splitlines():
            line = line.strip()
            if line:
                last = line
        if last:
            try:
                out[c] = int(last)
            except ValueError:
                pass
    return out


def run_qc(vcf_path: str, name: str) -> dict:
    errors: list[str] = []
    warnings: list[str] = []
    checks: dict = {}

    # --- Hard: file present ---
    if not os.path.isfile(vcf_path):
        errors.append(f"file does not exist: {vcf_path}")
        return _result(name, vcf_path, errors, warnings, checks)
    if os.path.getsize(vcf_path) == 0:
        errors.append(f"file is empty: {vcf_path}")
        return _result(name, vcf_path, errors, warnings, checks)
    checks["size_bytes"] = os.path.getsize(vcf_path)

    # --- Hard: tabix index ---
    tbi = vcf_path + ".tbi"
    csi = vcf_path + ".csi"
    if not (os.path.isfile(tbi) or os.path.isfile(csi)):
        errors.append(f"missing tabix index: expected {tbi} (or .csi)")
    checks["tabix_index_present"] = os.path.isfile(tbi) or os.path.isfile(csi)

    # --- Hard: header parseable ---
    header_proc = _run(["bcftools", "view", "-h", vcf_path])
    if header_proc.returncode != 0:
        errors.append(
            f"bcftools cannot read header: {header_proc.stderr.strip()[:200]}"
        )
        return _result(name, vcf_path, errors, warnings, checks)
    parsed = _parse_header(header_proc.stdout)
    checks["header_parseable"] = True
    checks["reference"] = parsed["reference"]
    checks["info_fields"] = parsed["info_ids"]
    checks["has_format_gt"] = parsed["has_format_gt"]
    checks["sample_count"] = len(parsed["samples"])

    # --- Hard: #CHROM line + samples ---
    if not parsed["chrom_line_present"]:
        errors.append("header has no #CHROM line")
    elif not parsed["samples"]:
        errors.append("#CHROM line has no sample columns")

    # --- Hard: indexed contigs + record count ---
    idx_proc = _run(["bcftools", "index", "-s", vcf_path])
    contigs, total_records = _parse_index(idx_proc.stdout)
    checks["contigs"] = contigs
    checks["n_variants"] = total_records
    if not contigs:
        errors.append("index reports no contigs (file may have no data)")
    if total_records is not None and total_records == 0:
        errors.append("file contains zero variant records")
    elif total_records is None and contigs:
        # Older indices lack counts — fall back to a streaming peek for the
        # presence of at least one record.
        r = _run(["bcftools", "view", "-H", "-r", contigs[0], vcf_path])
        if not r.stdout.strip():
            errors.append("file contains no variant records in first contig")

    # --- Soft: reference build ---
    build = _detect_build(parsed["reference"])
    checks["reference_build"] = build
    if not build:
        warnings.append(
            f"reference '{parsed['reference'] or '(none)'}' does not match "
            "any recognised human build (hg19/GRCh37/hg38/GRCh38)"
        )

    # --- Soft: human contigs ---
    human_contigs = [c for c in contigs if c in HUMAN_CHROMS]
    non_human = [c for c in contigs if c not in HUMAN_CHROMS]
    checks["human_contigs"] = human_contigs
    checks["non_human_contigs"] = non_human
    if contigs and not human_contigs:
        warnings.append(
            f"no indexed contigs match human chromosome names; saw {contigs}"
        )
    elif non_human:
        warnings.append(
            f"indexed contigs include names that are not standard human "
            f"chromosomes: {non_human}"
        )

    # --- Soft: FORMAT=GT ---
    if not parsed["has_format_gt"]:
        warnings.append("header does not declare FORMAT=GT")

    # --- Soft: INFO fields the downstream pipeline relies on ---
    has_af = "AF" in parsed["info_ids"]
    has_ac_an = "AC" in parsed["info_ids"] and "AN" in parsed["info_ids"]
    checks["info_has_af"] = has_af
    checks["info_has_ac_an"] = has_ac_an
    if not has_af and not has_ac_an:
        warnings.append(
            "header declares neither AF nor AC+AN — scan_variant will not be "
            "able to classify allele frequency for this file"
        )

    # --- Soft: variant position sanity ---
    if contigs and not errors:
        # Only when the rest of the file checks out — otherwise don't spend
        # the I/O. One-contig-per-file (1000G layout) keeps this cheap.
        max_pos = _max_position_per_contig(vcf_path, contigs)
        checks["max_position_per_contig"] = max_pos
        overlong = {c: p for c, p in max_pos.items()
                    if p > MAX_HUMAN_CHROM_LENGTH}
        if overlong:
            warnings.append(
                f"variant positions exceed plausible human chromosome "
                f"length ({MAX_HUMAN_CHROM_LENGTH:,} bp): {overlong}"
            )

    return _result(name, vcf_path, errors, warnings, checks)


def _result(name: str, vcf_path: str, errors: list[str],
            warnings: list[str], checks: dict) -> dict:
    return {
        "name": name,
        "file": os.path.basename(vcf_path),
        "pass": not errors,
        "errors": errors,
        "warnings": warnings,
        "checks": checks,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--vcf", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--out", required=True, help="Output JSON path")
    p.add_argument("--strict", action="store_true",
                   help="Exit 1 if any hard check fails")
    args = p.parse_args()

    result = run_qc(args.vcf, args.name)

    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)

    if result["errors"]:
        print(f"[qc_validate] HARD FAILURES for {args.name}:", file=sys.stderr)
        for e in result["errors"]:
            print(f"  - {e}", file=sys.stderr)
    if result["warnings"]:
        print(f"[qc_validate] WARNINGS for {args.name}:", file=sys.stderr)
        for w in result["warnings"]:
            print(f"  - {w}", file=sys.stderr)

    if args.strict and result["errors"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
