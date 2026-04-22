#!/usr/bin/env python3
"""Look for a variant (chrom+pos+REF+ALT) in a VCF and emit per-sample carriers.

Produces two files named from --out-prefix:
  <prefix>.variant.json   — scan result (status, AF, counts)
  <prefix>.carriers.tsv   — one row per carrier sample; header-only if none

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


CARRIERS_HEADER = "file\tsample\tvariant_id\tchrom\tpos\tref\talt\tgenotype\talt_dosage"


def run(cmd, shell=False):
    return subprocess.run(cmd, capture_output=True, text=True,
                          check=False, shell=shell).stdout


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def dosage(gt):
    return sum(1 for c in gt if c == "1")


def write_carriers(path, rows):
    with open(path, "w") as f:
        f.write(CARRIERS_HEADER + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--pos", type=int, required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--alt", required=True)
    p.add_argument("--min-af", type=float, default=0.0)
    p.add_argument("--max-af", type=float, default=1.0)
    args = p.parse_args()

    json_path = f"{args.out_prefix}.variant.json"
    carriers_path = f"{args.out_prefix}.carriers.tsv"
    file_base = os.path.basename(args.vcf)

    result = {
        "name": args.name,
        "file": file_base,
        "variant_chrom": args.chrom,
        "variant_pos": args.pos,
        "variant_ref": args.ref,
        "variant_alt": args.alt,
        "min_af": args.min_af,
        "max_af": args.max_af,
    }

    def finish(status, **extra):
        result["status"] = status
        result.update(extra)
        with open(json_path, "w") as f:
            json.dump(result, f, indent=2)

    idx = run(["bcftools", "index", "-s", args.vcf])
    contigs = [l.split()[0] for l in idx.strip().splitlines() if l.strip()]
    if not contigs:
        contigs = [c for c in run(["tabix", "-l", args.vcf]).splitlines() if c]

    bare = args.chrom.replace("chr", "")
    prefixed = f"chr{bare}"
    matched_contig = next((c for c in (args.chrom, bare, prefixed) if c in contigs), None)

    if not matched_contig:
        write_carriers(carriers_path, [])
        finish("not_applicable",
               note=f"chromosome '{args.chrom}' not present in VCF")
        return

    region = f"{matched_contig}:{args.pos}-{args.pos}"
    view = run(["bcftools", "view", "-r", region, "-H", args.vcf])
    rows = [r.split("\t") for r in view.strip().splitlines() if r.strip()]
    matches = [r for r in rows if r[3] == args.ref and r[4] == args.alt]

    if not matches:
        write_carriers(carriers_path, [])
        if not rows:
            finish("position_empty",
                   note=f"no variant calls at {matched_contig}:{args.pos}")
        else:
            seen = ",".join(f"{r[3]}>{r[4]}" for r in rows)
            finish("absent",
                   note=f"no {args.ref}>{args.alt} at position (saw: {seen})")
        return

    record = matches[0]
    info = {}
    for kv in record[7].split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            info[k] = v

    af = to_float(info.get("AF"))
    if af is None:
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

    # Per-sample genotypes at this site
    gt_out = run([
        "bcftools", "query",
        "-r", region,
        "-i", f'REF="{args.ref}" & ALT="{args.alt}"',
        "-f", "[%SAMPLE\\t%GT\\n]",
        args.vcf,
    ])

    variant_id = record[2]
    chrom = record[0]
    pos = record[1]
    ref = record[3]
    alt = record[4]

    carriers = []
    n_samples_scanned = 0
    het = hom = 0
    for line in gt_out.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        sample, gt = parts
        n_samples_scanned += 1
        d = dosage(gt)
        if d >= 1:
            carriers.append((file_base, sample, variant_id, chrom, pos,
                             ref, alt, gt, d))
            if d == 1:
                het += 1
            elif d == 2:
                hom += 1

    write_carriers(carriers_path, carriers)

    finish(status,
           variant_id=variant_id,
           chrom=chrom, pos=int(pos), ref=ref, alt=alt,
           ac=to_float(info.get("AC")),
           an=to_float(info.get("AN")),
           af=af,
           eas_af=to_float(info.get("EAS_AF")),
           eur_af=to_float(info.get("EUR_AF")),
           afr_af=to_float(info.get("AFR_AF")),
           sas_af=to_float(info.get("SAS_AF")),
           amr_af=to_float(info.get("AMR_AF")),
           n_samples_scanned=n_samples_scanned,
           n_carriers=len(carriers),
           n_heterozygotes=het,
           n_homozygous_alt=hom)


if __name__ == "__main__":
    main()
