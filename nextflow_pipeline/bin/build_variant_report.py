#!/usr/bin/env python3
"""Aggregate per-VCF variant-scan JSONs into a markdown report."""

import argparse
import glob
import json
import os
from collections import Counter


def fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--variant-name", required=True)
    p.add_argument("--variant-chrom", required=True)
    p.add_argument("--variant-pos", type=int, required=True)
    p.add_argument("--variant-ref", required=True)
    p.add_argument("--variant-alt", required=True)
    p.add_argument("--min-af", type=float, required=True)
    p.add_argument("--max-af", type=float, required=True)
    args = p.parse_args()

    results = []
    for f in sorted(glob.glob(os.path.join(args.input_dir, "*.variant.json"))):
        with open(f) as fh:
            results.append(json.load(fh))

    status_counts = Counter(r.get("status", "unknown") for r in results)
    in_range = [r for r in results if r.get("status") == "present_in_range"]

    out = [f"# Variant Scan Report — {args.variant_name}", ""]
    out += ["## Variant specification", ""]
    out += [f"- **Name:** {args.variant_name}"]
    out += [f"- **Coordinates:** `{args.variant_chrom}:{args.variant_pos}`"]
    out += [f"- **Alleles:** `{args.variant_ref}` → `{args.variant_alt}`"]
    out += [f"- **Acceptable AF range:** [{args.min_af}, {args.max_af}]"]
    out += [""]

    out += ["## Summary", ""]
    out += [f"- **Files scanned:** {len(results)}"]
    for status, count in status_counts.most_common():
        out += [f"- **{status}:** {count}"]
    out += [""]

    out += [f"## Files with variant present within AF range ({len(in_range)})", ""]
    if in_range:
        out += ["| File | Variant ID | REF | ALT | AF | EUR_AF | AFR_AF | EAS_AF | SAS_AF | AMR_AF |"]
        out += ["|---|---|---|---|---|---|---|---|---|---|"]
        for r in in_range:
            out += [
                f"| `{r['file']}` | {r.get('variant_id', '.')} "
                f"| {r.get('ref', '-')} | {r.get('alt', '-')} "
                f"| {fmt(r.get('af'))} | {fmt(r.get('eur_af'))} | {fmt(r.get('afr_af'))} "
                f"| {fmt(r.get('eas_af'))} | {fmt(r.get('sas_af'))} | {fmt(r.get('amr_af'))} |"
            ]
    else:
        out += ["_No files contain the variant within the specified AF range._"]
    out += [""]

    out += ["## Full results", ""]
    out += ["| File | Status | Note | AF |"]
    out += ["|---|---|---|---|"]
    for r in results:
        out += [
            f"| `{r['file']}` | {r.get('status', '-')} "
            f"| {r.get('note', '')} | {fmt(r.get('af'))} |"
        ]

    with open(args.output, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
