"""Write a single-sample bgzipped + tabix-indexed VCF."""

from __future__ import annotations

import random
import subprocess
from pathlib import Path

from .background import alt_dosages
from .builds import BUILDS
from .errors import (
    DEFAULT_DROPOUT_RATE,
    DEFAULT_GT_ERROR_RATE,
    maybe_dropout,
    maybe_flip_gt,
    new_error_stats,
)
from .header import build_header
from .quality import (
    DEFAULT_DP_MEAN,
    DEFAULT_DP_SAMPLE_JITTER_SD,
    draw_site_quality,
    gq_from_ad,
    sample_lambda,
)


def _contig_sort_key(chrom: str, pos: int,
                     contig_order: dict) -> tuple:
    """Order chromosomes using the reference dict; unknowns sort last."""
    return (contig_order.get(chrom, len(contig_order)), pos)


def write_person_vcf(out_path: Path, person: dict, build: str,
                     rng: random.Random,
                     dp_mean: float = DEFAULT_DP_MEAN,
                     error_rate: float = 0.0,
                     dropout_rate: float = 0.0,
                     stats: dict | None = None) -> Path:
    """Write a single-sample bgzipped+indexed VCF.

    Each record carries simulated DP/GQ/AD alongside GT. The per-sample
    mean depth is drawn once (Gaussian around `dp_mean`) so different
    samples in a cohort show plausibly different coverage profiles.

    M9 sequencing-error parameters:
      * `error_rate` — per-call probability of a GT flip applied
        *after* AD has been drawn from the truth, so the recomputed GQ
        naturally drops where the call disagrees with the reads.
      * `dropout_rate` — per-call probability of a coverage dropout
        (DP=0, AD all-zero, GQ=0, GT=`./.`).
      * `stats` — optional dict mutated in place with running counts of
        `flipped`, `dropped`, `total_calls`. The caller seeds it with
        `errors.new_error_stats()`.
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
            alts = variant["alts"]
            n_alts = len(alts)
            n_alleles = n_alts + 1  # +1 for REF

            # AC is per-alt (Number=A); AN is total allele count (fixed at 2
            # for a diploid single-sample). AF = AC / AN per alt.
            per_alt_dosage = alt_dosages(gt, n_alts)
            ac_str = ",".join(str(d) for d in per_alt_dosage)
            af_str = ",".join(f"{d/2:.1f}" for d in per_alt_dosage)
            info_parts = [f"AC={ac_str}", "AN=2", f"AF={af_str}"]
            if is_hi:
                info_parts.append("HIGHLIGHT")
            # Annotation overlays (M7) are picked up regardless of
            # whether this is the highlighted record. Cohort-level
            # ClinVar / COSMIC injection lands the same INFO tags on
            # background records, so the writer carries them over here.
            if variant.get("clnsig") and variant["clnsig"] != ".":
                info_parts.append(f"CLNSIG={variant['clnsig']}")
            if variant.get("clndn") and variant["clndn"] != ".":
                info_parts.append(f"CLNDN={variant['clndn']}")
            if variant.get("cosmic_id"):
                info_parts.append(f"COSMIC_ID={variant['cosmic_id']}")
            if variant.get("cosmic_gene"):
                info_parts.append(f"COSMIC_GENE={variant['cosmic_gene']}")
            # Structural variants (M8) declare SVTYPE / SVLEN / END /
            # CIPOS. The ALT is a symbolic <DEL>/<DUP>/<INV> already in
            # variant["alts"], so the standard ",".join(alts) writes
            # the right thing — only INFO needs extra fields.
            if variant.get("svtype"):
                info_parts.append(f"SVTYPE={variant['svtype']}")
                if variant.get("svlen") is not None:
                    info_parts.append(f"SVLEN={variant['svlen']}")
                if variant.get("end") is not None:
                    info_parts.append(f"END={variant['end']}")
                if variant.get("cipos"):
                    lo, hi = variant["cipos"]
                    info_parts.append(f"CIPOS={lo},{hi}")

            # Truth-state DP/AD/GQ first; M9 noise is layered on top so
            # AD reflects what the reads would have looked like under
            # the *true* genotype and any GT flip has to live with low
            # GQ (because gq_from_ad sees the inconsistency).
            dp, ad, gq = draw_site_quality(gt, n_alleles, sample_lam, rng)
            called_gt = gt
            if stats is not None:
                stats["total_calls"] = stats.get("total_calls", 0) + 1
            if dropout_rate > 0 and maybe_dropout(rng, dropout_rate):
                called_gt = "./."
                dp = 0
                ad = (0,) * n_alleles
                gq = 0
                if stats is not None:
                    stats["dropped"] = stats.get("dropped", 0) + 1
            elif error_rate > 0:
                new_gt, flipped = maybe_flip_gt(gt, rng, error_rate)
                if flipped:
                    called_gt = new_gt
                    gq = gq_from_ad(called_gt, ad)
                    if stats is not None:
                        stats["flipped"] = stats.get("flipped", 0) + 1
            ad_str = ",".join(str(x) for x in ad)
            sample_field = f"{called_gt}:{dp}:{gq}:{ad_str}"
            fh.write("\t".join([
                variant["chrom"], str(variant["pos"]),
                variant.get("id") or ".",
                variant["ref"], ",".join(alts),
                "100", "PASS", ";".join(info_parts),
                "GT:DP:GQ:AD", sample_field,
            ]) + "\n")

    subprocess.run(["bgzip", "-f", str(plain_path)],
                   check=True, capture_output=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", str_out],
                   check=True, capture_output=True)
    return out_path
