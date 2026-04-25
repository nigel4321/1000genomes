"""Sequencing-error injection (M9).

Lightweight model that perturbs per-call genotypes and read coverage at
a configured target false-discovery rate. The perturbation is applied
**after** AD/DP have been drawn from the true genotype, so the called
GT can disagree with what the reads actually support. This produces the
realistic behaviour that GQ drops when the caller "gets it wrong" — the
recomputed GQ from gq_from_ad sees that AD is inconsistent with the
emitted GT and reports low confidence.

Two failure modes are modelled:

* **GT flips** — at probability `error_rate`, the called GT is swapped
  for a different value. For biallelic calls (`0|0`, `0|1`, `1|0`,
  `1|1`) the flip distribution mirrors empirical caller behaviour: hom
  → het is more common than hom → opposite-hom; het splits between
  hom-ref and hom-alt. Multi-allelic hets (e.g. `1|2`) collapse one
  allele to REF.
* **Dropouts** — at probability `dropout_rate`, the entire call drops
  coverage (DP=0, AD=zeros, GQ=0, GT=`./.`). Models stochastic regions
  of low/no read support.

Heavy-path ART read simulation is gated behind a separate `--art` flag
and not implemented in this milestone — that path will land alongside
the M11 GRCh38 reference FASTA.
"""

from __future__ import annotations

import random


# ~0.1% FDR target — realistic for a well-tuned germline caller after
# standard hard filters. Light enough that downstream stats (Ti/Tv,
# Het/Hom) are barely perturbed but visible in a truth-set comparison.
DEFAULT_GT_ERROR_RATE = 0.001
# Dropouts are rarer than mis-calls in a well-covered modern WGS run;
# 0.05% is a conservative default that still surfaces a handful per
# 10k-record VCF so M11's truth-set BED has something to grade.
DEFAULT_DROPOUT_RATE = 0.0005

# Distribution of biallelic flips. Empirically, callers under-call hets
# more often than they fabricate them; hom→het is the most common
# specific direction (low coverage / contamination), so we weight it
# higher than hom→opposite-hom.
_BIALLELIC_FLIPS = {
    "0|0": (("0|1", 0.7), ("1|1", 0.3)),
    "0|1": (("0|0", 0.5), ("1|1", 0.5)),
    "1|0": (("0|0", 0.5), ("1|1", 0.5)),
    "1|1": (("0|1", 0.7), ("0|0", 0.3)),
}


def _weighted_pick(options: tuple, rng: random.Random) -> str:
    r = rng.random()
    cum = 0.0
    for choice, w in options:
        cum += w
        if r < cum:
            return choice
    return options[-1][0]


def maybe_flip_gt(true_gt: str, rng: random.Random,
                  error_rate: float) -> tuple[str, bool]:
    """Return `(called_gt, flipped)`.

    With probability `error_rate`, replaces the called GT with a
    plausible mis-call. For biallelic GTs the swap follows
    `_BIALLELIC_FLIPS`; for `1|2`-style multi-allelic hets we collapse
    one allele to REF (the dominant real-world failure mode for
    multi-allelic sites). Anything we can't parse passes through
    unchanged.

    `error_rate <= 0` is a fast no-op.
    """
    if error_rate <= 0:
        return true_gt, False
    if rng.random() >= error_rate:
        return true_gt, False

    if true_gt in _BIALLELIC_FLIPS:
        return _weighted_pick(_BIALLELIC_FLIPS[true_gt], rng), True

    parts = true_gt.split("|")
    if len(parts) != 2:
        return true_gt, False  # haploid / phased-malformed → leave alone
    a, b = parts
    try:
        ai, bi = int(a), int(b)
    except ValueError:
        return true_gt, False
    if ai == bi:
        # Hom on a non-REF non-ALT-1 allele (rare): drop one allele.
        return f"0|{b}", True
    # Multi-allelic het like "1|2" → caller usually drops to "0|<one of them>"
    if rng.random() < 0.5:
        return f"0|{b}", True
    return f"{a}|0", True


def maybe_dropout(rng: random.Random, dropout_rate: float) -> bool:
    """`True` if this call should lose coverage entirely."""
    if dropout_rate <= 0:
        return False
    return rng.random() < dropout_rate


def new_error_stats() -> dict:
    """Fresh accumulator for per-batch error counters."""
    return {"flipped": 0, "dropped": 0, "total_calls": 0}


def merge_stats(into: dict, other: dict) -> None:
    """Combine `other` into `into` in place. Keys missing from `into`
    are seeded so partial accumulators add cleanly."""
    for k, v in other.items():
        into[k] = into.get(k, 0) + v
