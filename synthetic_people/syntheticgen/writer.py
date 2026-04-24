"""Write a single-sample bgzipped + tabix-indexed VCF."""

from __future__ import annotations

import random
import subprocess
from pathlib import Path

from .background import alt_dosage
from .builds import BUILDS
from .header import build_header
from .quality import (
    DEFAULT_DP_MEAN,
    DEFAULT_DP_SAMPLE_JITTER_SD,
    draw_site_quality,
    sample_lambda,
)


def _contig_sort_key(chrom: str, pos: int,
                     contig_order: dict) -> tuple:
    """Order chromosomes using the reference dict; unknowns sort last."""
    return (contig_order.get(chrom, len(contig_order)), pos)


def write_person_vcf(out_path: Path, person: dict, build: str,
                     rng: random.Random,
                     dp_mean: float = DEFAULT_DP_MEAN) -> Path:
    """Write a single-sample bgzipped+indexed VCF.

    Each record carries simulated DP/GQ/AD alongside GT. The per-sample
    mean depth is drawn once (Gaussian around `dp_mean`) so different
    samples in a cohort show plausibly different coverage profiles.
    """
    contigs = BUILDS[build]["contigs"]
    contig_order = {c: i for i, c in enumerate(contigs)}

    records: list[tuple[dict, str, bool]] = []
    records.append((person["highlighted"], person["highlighted"]["gt"], True))
    for bg in person["background"]:
        records.append((bg, bg["gt"], False))
    records.sort(key=lambda r: _contig_sort_key(
        r[0]["chrom"], r[0]["pos"], contig_order))

    str_out = str(out_path)
    if not str_out.endswith(".vcf.gz"):
        raise ValueError("out_path must end in .vcf.gz")
    plain_path = Path(str_out[:-len(".gz")])  # drop .gz → .vcf

    sample_lam = sample_lambda(dp_mean, DEFAULT_DP_SAMPLE_JITTER_SD, rng)

    with open(plain_path, "w") as fh:
        fh.write(build_header(build, person["sample_id"]))
        for variant, gt, is_hi in records:
            dosage = alt_dosage(gt)
            info_parts = [f"AC={dosage}", "AN=2", f"AF={dosage/2:.1f}"]
            if is_hi:
                info_parts.append("HIGHLIGHT")
                if variant.get("clnsig") and variant["clnsig"] != ".":
                    info_parts.append(f"CLNSIG={variant['clnsig']}")
                if variant.get("clndn") and variant["clndn"] != ".":
                    info_parts.append(f"CLNDN={variant['clndn']}")
            dp, ref_d, alt_d, gq = draw_site_quality(gt, sample_lam, rng)
            sample_field = f"{gt}:{dp}:{gq}:{ref_d},{alt_d}"
            fh.write("\t".join([
                variant["chrom"], str(variant["pos"]),
                variant.get("id") or ".",
                variant["ref"], variant["alt"],
                "100", "PASS", ";".join(info_parts),
                "GT:DP:GQ:AD", sample_field,
            ]) + "\n")

    subprocess.run(["bgzip", "-f", str(plain_path)],
                   check=True, capture_output=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", str_out],
                   check=True, capture_output=True)
    return out_path
