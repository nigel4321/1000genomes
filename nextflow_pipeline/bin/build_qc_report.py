#!/usr/bin/env python3
"""Aggregate per-VCF qc.json files into a markdown QC report."""

import argparse
import glob
import json
import os


def _status(pass_: bool) -> str:
    return "PASS" if pass_ else "FAIL"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    results = []
    for f in sorted(glob.glob(os.path.join(args.input_dir, "*.qc.json"))):
        with open(f) as fh:
            results.append(json.load(fh))

    n_total = len(results)
    n_pass = sum(1 for r in results if r.get("pass"))
    n_fail = n_total - n_pass
    n_warn = sum(1 for r in results if r.get("warnings"))

    out = ["# VCF QC Report", ""]
    out += ["## Summary", ""]
    out += [f"- **Files scanned:** {n_total}"]
    out += [f"- **Passed:** {n_pass}"]
    out += [f"- **Failed:** {n_fail}"]
    out += [f"- **Passed with warnings:** {n_warn - n_fail if n_warn > n_fail else n_warn}"]
    out += [""]

    if n_fail:
        out += ["## Hard failures", ""]
        out += [
            "_These files could not be validated and would break downstream "
            "stages. The pipeline aborts in strict mode._",
            "",
        ]
        for r in results:
            if r.get("pass"):
                continue
            out += [f"### `{r['file']}`", ""]
            for e in r.get("errors", []):
                out += [f"- {e}"]
            out += [""]

    warned = [r for r in results if r.get("warnings")]
    if warned:
        out += ["## Warnings", ""]
        out += [
            "_The file is usable but one or more attributes do not match "
            "typical human genomic VCFs — eyeball before trusting results._",
            "",
        ]
        for r in warned:
            out += [f"### `{r['file']}`", ""]
            for w in r["warnings"]:
                out += [f"- {w}"]
            out += [""]

    out += ["## Per-file detail", ""]
    out += ["| File | Status | Samples | Variants | Reference | Build | "
            "Contigs | AF/AC+AN | GT | Errors | Warnings |"]
    out += ["|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        c = r.get("checks", {})
        af_flag = (
            "AF" if c.get("info_has_af")
            else ("AC+AN" if c.get("info_has_ac_an") else "-")
        )
        gt_flag = "yes" if c.get("has_format_gt") else "no"
        contigs = c.get("contigs") or []
        contigs_str = ",".join(contigs) if contigs else "-"
        ref = c.get("reference") or "-"
        # Long reference URLs make the table unreadable — trim them.
        if len(ref) > 40:
            ref = ref[:37] + "..."
        out += [
            f"| `{r['file']}` | {_status(r.get('pass', False))} "
            f"| {c.get('sample_count', '-')} "
            f"| {c.get('n_variants', '?') if c.get('n_variants') is not None else '?'} "
            f"| {ref} "
            f"| {c.get('reference_build') or '-'} "
            f"| {contigs_str} "
            f"| {af_flag} "
            f"| {gt_flag} "
            f"| {len(r.get('errors', []))} "
            f"| {len(r.get('warnings', []))} |"
        ]

    with open(args.output, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
