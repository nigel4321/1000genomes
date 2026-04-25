"""Structural-variant generator (M8).

Produces a handful of per-person SVs (deletions, tandem duplications,
inversions) following the VCF 4.2 symbolic-allele convention:

* `REF` is a single anchor base at the SV's `POS` (we synthesise a
  random base — the real reference isn't on disk).
* `ALT` is a symbolic allele: `<DEL>`, `<DUP>`, or `<INV>`.
* `INFO`: `SVTYPE`, `SVLEN` (negative for DEL, positive for DUP/INV),
  `END` (POS + |SVLEN|), and `CIPOS` (a default ±50 bp confidence
  interval — every SV is "imprecise" at this stage).

Genotypes default to phased het with a small homozygous-alt tail. SV
loci are randomly placed across the simulated chromosome region,
biased away from `POS=1` so the anchor base + downstream span stays
in-bounds for the contig length recorded in the header.

The output dict shape mirrors the SNV `variant` dict consumed by
`writer.py` so SV records flow through the same write path with no
extra branching except for the symbolic-ALT INFO emission.
"""

from __future__ import annotations

import math
import random


SV_TYPES = ("DEL", "DUP", "INV")
DEFAULT_SV_TYPE_WEIGHTS = (0.50, 0.30, 0.20)

DEFAULT_SVS_PER_PERSON = 3
DEFAULT_SV_LENGTH_MIN_BP = 50
DEFAULT_SV_LENGTH_MAX_BP = 10_000

# Confidence-interval declared on every SV — we mark all of them
# imprecise because the simulated POS / END are random. Real callers
# emit tighter CIs for split-read-confirmed events; that's a future
# refinement, not in scope here.
DEFAULT_CIPOS_HALFWIDTH = 50


def _draw_length(rng: random.Random, lo_bp: int, hi_bp: int) -> int:
    """Log-uniform draw — SV length distributions are heavy-tailed."""
    if lo_bp < 1 or hi_bp < lo_bp:
        raise ValueError(f"bad SV length range: [{lo_bp}, {hi_bp}]")
    if lo_bp == hi_bp:
        return lo_bp
    log_lo = math.log(lo_bp)
    log_hi = math.log(hi_bp)
    return max(lo_bp, min(hi_bp, int(round(math.exp(
        rng.uniform(log_lo, log_hi))))))


def _draw_anchor_base(rng: random.Random) -> str:
    return rng.choice(("A", "C", "G", "T"))


def _draw_gt(rng: random.Random) -> str:
    """Phased SV genotype: ~80% het, ~20% hom-alt."""
    return rng.choices(("0|1", "1|1"), weights=(0.80, 0.20))[0]


def _build_sv_record(chrom: str, pos: int, length_bp: int,
                     svtype: str, gt: str,
                     rng: random.Random,
                     cipos_halfwidth: int = DEFAULT_CIPOS_HALFWIDTH
                     ) -> dict:
    """Assemble one SV variant dict in the writer-compatible shape."""
    if svtype not in SV_TYPES:
        raise ValueError(f"unknown SVTYPE {svtype!r}")
    end = pos + length_bp
    svlen = -length_bp if svtype == "DEL" else length_bp
    return {
        "chrom": chrom,
        "pos": pos,
        "id": ".",
        "ref": _draw_anchor_base(rng),
        "alts": [f"<{svtype}>"],
        "afs": [None],          # writer fills AF from the GT dosage
        "gt": gt,
        "svtype": svtype,
        "svlen": svlen,
        "end": end,
        "cipos": (-cipos_halfwidth, cipos_halfwidth),
    }


def generate_person_svs(rng: random.Random,
                        chromosomes: list[str],
                        chrom_length_bp: int,
                        n_svs: int = DEFAULT_SVS_PER_PERSON,
                        length_min_bp: int = DEFAULT_SV_LENGTH_MIN_BP,
                        length_max_bp: int = DEFAULT_SV_LENGTH_MAX_BP,
                        type_weights: tuple = DEFAULT_SV_TYPE_WEIGHTS,
                        ) -> list[dict]:
    """Produce `n_svs` SVs distributed uniformly across `chromosomes`.

    `chrom_length_bp` is the simulated span — SV positions are drawn
    from `[1, chrom_length_bp - length_max_bp]` so the entire SV
    (anchor base + downstream) stays within the simulated region. The
    chromosome list is sampled with replacement; tiny `n_svs` won't
    spread evenly across a large chrom set, which is fine.
    """
    if n_svs <= 0:
        return []
    if not chromosomes:
        raise ValueError("at least one chromosome required")
    if chrom_length_bp <= length_max_bp + 1:
        raise ValueError(
            f"chromosome span {chrom_length_bp} too small for SV "
            f"length up to {length_max_bp}"
        )

    pos_max = chrom_length_bp - length_max_bp
    out: list[dict] = []
    for _ in range(n_svs):
        chrom = rng.choice(chromosomes)
        length_bp = _draw_length(rng, length_min_bp, length_max_bp)
        pos = rng.randint(1, pos_max)
        svtype = rng.choices(SV_TYPES, weights=type_weights)[0]
        gt = _draw_gt(rng)
        out.append(_build_sv_record(chrom, pos, length_bp, svtype, gt,
                                    rng))
    return out
