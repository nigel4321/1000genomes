#!/usr/bin/env python3
"""Look for a specific variant (by chrom+pos+REF+ALT) in a VCF. Emits JSON.

Statuses:
  not_applicable          — target chromosome is not in this VCF
  position_empty          — chromosome is present but no variant call at position
  absent                  — variant(s) at position but none with matching REF/ALT
  present_in_range        — variant matches and AF in [min_af, max_af]
  present_below_threshold — variant matches but AF < min_af
  present_above_threshold — variant matches but AF > max_af
  present_af_unknown      — variant matches but AF could not be determined
"""

import argparse
import json
import os
import subprocess
import sys


def run(cmd, shell=False):
    return subprocess.run(cmd, capture_output=True, text=True, check=False, shell=shell).stdout


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--pos", type=int, required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--alt", required=True)
    p.add_argument("--min-af", type=float, default=0.0)
    p.add_argument("--max-af", type=float, default=1.0)
    args = p.parse_args()

    result = {
        "name": args.name,
        "file": os.path.basename(args.vcf),
        "variant_chrom": args.chrom,
        "variant_pos": args.pos,
        "variant_ref": args.ref,
        "variant_alt": args.alt,
        "min_af": args.min_af,
        "max_af": args.max_af,
    }

    idx = run(["bcftools", "index", "-s", args.vcf])
    contigs = [l.split()[0] for l in idx.strip().splitlines() if l.strip()]
    if not contigs:
        contigs = [c for c in run(["tabix", "-l", args.vcf]).splitlines() if c]

    bare = args.chrom.replace("chr", "")
    prefixed = f"chr{bare}"
    matched_contig = next((c for c in (args.chrom, bare, prefixed) if c in contigs), None)

    if not matched_contig:
        result.update(status="not_applicable",
                      note=f"chromosome '{args.chrom}' not present in VCF")
        print(json.dumps(result, indent=2))
        return

    region = f"{matched_contig}:{args.pos}-{args.pos}"
    view = run(["bcftools", "view", "-r", region, "-H", args.vcf])
    rows = [r.split("\t") for r in view.strip().splitlines() if r.strip()]
    matches = [r for r in rows if r[3] == args.ref and r[4] == args.alt]

    if not matches:
        if not rows:
            result.update(status="position_empty",
                          note=f"no variant calls at {matched_contig}:{args.pos}")
        else:
            seen = ",".join(f"{r[3]}>{r[4]}" for r in rows)
            result.update(status="absent",
                          note=f"no {args.ref}>{args.alt} at position (saw: {seen})")
        print(json.dumps(result, indent=2))
        return

    record = matches[0]
    info = {}
    for kv in record[7].split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            info[k] = v

    af = to_float(info.get("AF"))
    if af is None:
        # Recompute AF from genotypes if not pre-computed in INFO.
        cmd = (
            f"bcftools +fill-tags {args.vcf} -r {region} -- -t AF | "
            f"bcftools query -i 'REF=\"{args.ref}\" & ALT=\"{args.alt}\"' "
            f"-f '%INFO/AF\\n'"
        )
        filled = run(cmd, shell=True).strip().splitlines()
        af = to_float(filled[0]) if filled else None

    if af is None:
        status = "present_af_unknown"
    elif af < args.min_af:
        status = "present_below_threshold"
    elif af > args.max_af:
        status = "present_above_threshold"
    else:
        status = "present_in_range"

    result.update({
        "status": status,
        "variant_id": record[2],
        "chrom": record[0],
        "pos": int(record[1]),
        "ref": record[3],
        "alt": record[4],
        "ac": to_float(info.get("AC")),
        "an": to_float(info.get("AN")),
        "af": af,
        "eas_af": to_float(info.get("EAS_AF")),
        "eur_af": to_float(info.get("EUR_AF")),
        "afr_af": to_float(info.get("AFR_AF")),
        "sas_af": to_float(info.get("SAS_AF")),
        "amr_af": to_float(info.get("AMR_AF")),
    })
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
