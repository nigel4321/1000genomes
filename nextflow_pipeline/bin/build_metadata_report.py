#!/usr/bin/env python3
"""Aggregate per-VCF metadata JSONs into a markdown cohort overview."""

import argparse
import glob
import json
import os
from collections import Counter


def human_size(n):
    f = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if f < 1024:
            return f"{f:.1f} {unit}"
        f /= 1024
    return f"{f:.1f} PB"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    metas = []
    for f in sorted(glob.glob(os.path.join(args.input_dir, "*.meta.json"))):
        with open(f) as fh:
            metas.append(json.load(fh))

    all_contigs = sorted({c for m in metas for c in m.get("contigs", [])},
                         key=lambda x: (len(x), x))
    total_variants = sum(m.get("n_variants") or 0 for m in metas)
    any_unknown_counts = any(m.get("n_variants") is None for m in metas)
    total_size = sum(m.get("size_bytes", 0) for m in metas)
    cohorts = Counter(m.get("sample_hash") for m in metas)
    references = Counter(m.get("reference") or "unknown" for m in metas)
    pipelines = Counter(m.get("pipeline_tag") or "(none)" for m in metas)

    out = ["# VCF Cohort Metadata Report", ""]
    out += ["## Overview", ""]
    out += [f"- **Files scanned:** {len(metas)}"]
    variants_note = " (partial — some files had no count metadata in their index)" if any_unknown_counts else ""
    out += [f"- **Total variants across files:** {total_variants:,}{variants_note}"]
    out += [f"- **Total size on disk:** {human_size(total_size)}"]
    out += [f"- **Contigs covered ({len(all_contigs)}):** {', '.join(all_contigs) or '-'}"]
    out += [f"- **Distinct sample cohorts (by sample-list hash):** {len(cohorts)}"]
    if len(cohorts) == 1 and metas:
        out += [f"  - All files share the same {metas[0]['n_samples']:,}-sample cohort "
                "(likely a per-chromosome split of one dataset)."]
    elif len(cohorts) > 1:
        out += ["  - Sample sets differ between files — downstream joins must account for this."]
    out += [f"- **Reference builds:** " +
            ", ".join(f"{k} ({v})" for k, v in references.most_common())]
    out += [f"- **Pipeline tags:** " +
            ", ".join(f"{k} ({v})" for k, v in pipelines.most_common())]
    out += [""]

    out += ["## Per-file details", ""]
    out += ["| File | Contigs | Samples | Variants | Size | Reference | Pipeline | Date |"]
    out += ["|---|---|---|---|---|---|---|---|"]
    for m in metas:
        contigs = ",".join(m.get("contigs", [])) or "-"
        nv = m.get("n_variants")
        nv_str = f"{nv:,}" if nv is not None else "?"
        out += [
            f"| `{m['file']}` "
            f"| {contigs} "
            f"| {m['n_samples']:,} "
            f"| {nv_str} "
            f"| {human_size(m['size_bytes'])} "
            f"| {m.get('reference') or '-'} "
            f"| {m.get('pipeline_tag') or '-'} "
            f"| {m.get('date_stamp') or '-'} |"
        ]

    with open(args.output, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
