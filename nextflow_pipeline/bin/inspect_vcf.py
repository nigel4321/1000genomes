#!/usr/bin/env python3
"""Collect lightweight metadata for a single VCF via bcftools. Emits JSON."""

import hashlib
import json
import os
import re
import subprocess
import sys


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <vcf> <name>")
    vcf, name = sys.argv[1], sys.argv[2]

    idx = run(["bcftools", "index", "-s", vcf])
    contigs, n_variants = [], None
    for line in idx.strip().splitlines():
        parts = line.split()
        if len(parts) >= 3:
            contigs.append(parts[0])
            try:
                n_variants = (n_variants or 0) + int(parts[2])
            except ValueError:
                pass
    # Older tabix indices lack count metadata — fall back to `tabix -l`
    # for contigs and leave variant count as null.
    if not contigs:
        contigs = [c for c in run(["tabix", "-l", vcf]).splitlines() if c]

    samples = [s for s in run(["bcftools", "query", "-l", vcf]).splitlines() if s]
    sample_hash = hashlib.md5("\n".join(sorted(samples)).encode()).hexdigest()

    reference = "unknown"
    for line in run(["bcftools", "view", "-h", vcf]).splitlines():
        if line.startswith("##reference"):
            reference = line.split("=", 1)[1].strip()
            break

    base = os.path.basename(vcf)
    pipeline_match = re.search(r"phase\d+_[a-z0-9_]+", base)
    date_match = re.search(r"\d{8}", base)

    print(json.dumps({
        "name": name,
        "file": base,
        "contigs": contigs,
        "n_variants": n_variants,
        "n_samples": len(samples),
        "sample_hash": sample_hash,
        "reference": reference,
        "size_bytes": os.path.getsize(vcf),
        "pipeline_tag": pipeline_match.group(0) if pipeline_match else "",
        "date_stamp": date_match.group(0) if date_match else "",
    }, indent=2))


if __name__ == "__main__":
    main()
