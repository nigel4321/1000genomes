#!/usr/bin/env python3
"""Aggregate per-file carrier TSVs into a master TSV and markdown summary."""

import argparse
import glob
import os
from collections import defaultdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--output-md", required=True)
    p.add_argument("--variant-name", required=True)
    args = p.parse_args()

    header = None
    rows = []
    per_file = defaultdict(lambda: {"het": 0, "hom": 0})

    for f in sorted(glob.glob(os.path.join(args.input_dir, "*.carriers.tsv"))):
        with open(f) as fh:
            lines = fh.read().splitlines()
        if not lines:
            continue
        if header is None:
            header = lines[0]
        for line in lines[1:]:
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            rows.append(fields)
            try:
                d = int(fields[8])
                key = fields[0]
                if d == 1:
                    per_file[key]["het"] += 1
                elif d == 2:
                    per_file[key]["hom"] += 1
            except ValueError:
                pass

    with open(args.output_tsv, "w") as f:
        f.write((header or "file\tsample\tvariant_id\tchrom\tpos\tref\talt\tgenotype\talt_dosage") + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    total = len(rows)
    het = sum(1 for r in rows if r[8] == "1")
    hom = sum(1 for r in rows if r[8] == "2")
    unique_samples = len({r[1] for r in rows})
    files_with_data = len(per_file)

    out = [f"# Carrier Report — {args.variant_name}", ""]
    out += ["## Summary", ""]
    out += [f"- **Carrier records (rows in `carriers.tsv`):** {total}"]
    out += [f"- **Unique individuals carrying the alt allele:** {unique_samples}"]
    out += [f"- **Heterozygotes (dosage = 1):** {het}"]
    out += [f"- **Homozygous alt (dosage = 2):** {hom}"]
    out += [f"- **Files containing the variant:** {files_with_data}"]
    if total:
        ac_check = het + 2 * hom
        out += [f"- **Allele-count integrity check:** het + 2·hom = {het} + 2·{hom} = **{ac_check}**"]
    out += [""]

    if per_file:
        out += ["## Per-file breakdown", ""]
        out += ["| File | Heterozygotes | Homozygous alt | Total carriers |"]
        out += ["|---|---|---|---|"]
        for fname in sorted(per_file):
            c = per_file[fname]
            out += [f"| `{fname}` | {c['het']} | {c['hom']} | {c['het'] + c['hom']} |"]
        out += [""]

    out += ["## Output files", ""]
    out += ["- `carriers.tsv` — one row per (file, sample) where the sample carries the alt allele."]
    out += ["  - Columns: `file, sample, variant_id, chrom, pos, ref, alt, genotype, alt_dosage`"]
    out += [""]

    if total:
        out += ["## First 10 carriers (preview)", ""]
        out += ["| sample | genotype | alt_dosage |"]
        out += ["|---|---|---|"]
        for r in rows[:10]:
            out += [f"| {r[1]} | {r[7]} | {r[8]} |"]
        out += [""]

    with open(args.output_md, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
