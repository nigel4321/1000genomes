"""Site-frequency-spectrum sampler — draws minor allele counts ∝ 1/k^α.

The neutral, constant-size coalescent gives expected SFS proportional to
1/k (Watterson). Real human populations have undergone strong growth, so
empirical spectra are steeper than 1/k — singletons (AC=1) dominate. We
draw from `P(k) ∝ 1/k^α` with α configurable; the default α=2.0 reproduces
the gnomAD-style "more than half of variants are singletons" behaviour
across a wide range of cohort sizes.

This module underpins the cohort-level site generation in `cohort.py` and
is the source of truth for singleton fractions reported in the run log.
"""

from __future__ import annotations

import random
from pathlib import Path


DEFAULT_SFS_ALPHA = 2.0


def _weights(n_haplotypes: int, alpha: float) -> list[float]:
    return [1.0 / (k ** alpha) for k in range(1, n_haplotypes)]


def draw_minor_count(n_haplotypes: int, alpha: float,
                     rng: random.Random) -> int:
    """Draw a minor allele count k in {1, ..., n_haplotypes-1} with P(k) ∝ 1/k^α."""
    if n_haplotypes < 2:
        raise ValueError(f"n_haplotypes must be ≥ 2, got {n_haplotypes}")
    if alpha <= 0:
        raise ValueError(f"alpha must be positive, got {alpha}")
    ks = list(range(1, n_haplotypes))
    return rng.choices(ks, weights=_weights(n_haplotypes, alpha), k=1)[0]


def draw_allele_counts(n_haplotypes: int, n_alts: int, alpha: float,
                       rng: random.Random,
                       max_attempts: int = 50) -> list[int]:
    """Per-alt count vector drawn from the SFS with total ≤ n_haplotypes-1.

    Each alt is drawn independently from the 1/k^α distribution. If the
    total exceeds n_haplotypes-1 (which would require more alt alleles
    than haplotype slots while leaving at least one REF), we re-draw up
    to `max_attempts` times, then fall back to shrinking the largest
    counts by one until the total fits. The "≥ 1 REF" constraint keeps
    every cohort site truly variable (never fixed for alt).
    """
    if n_alts < 1:
        raise ValueError(f"n_alts must be ≥ 1, got {n_alts}")
    max_total = n_haplotypes - 1
    for _ in range(max_attempts):
        counts = [draw_minor_count(n_haplotypes, alpha, rng)
                  for _ in range(n_alts)]
        if sum(counts) <= max_total:
            return counts
    # Fallback: take the last draw and shrink the largest counts until
    # the total fits. This biases the result away from the raw SFS but
    # only kicks in for degenerate inputs (many alts + small cohort).
    total = sum(counts)
    while total > max_total:
        idx = max(range(n_alts), key=lambda k: counts[k])
        counts[idx] -= 1
        total -= 1
        if counts[idx] < 1:
            counts[idx] = 1
            # Shrinking below 1 means the inputs are unsatisfiable; bail
            # by proportional rescaling instead.
            scale = max_total / sum(counts)
            counts = [max(1, int(c * scale)) for c in counts]
            # Final trim.
            while sum(counts) > max_total:
                i = max(range(n_alts), key=lambda k: counts[k])
                if counts[i] > 1:
                    counts[i] -= 1
                else:
                    break
            break
    return counts


def sfs_histogram(sites: list) -> dict:
    """Return {k: n_alt_observations_with_count_k}.

    For multi-allelic sites, each alt contributes its own AC to the
    histogram independently — so a site with ALTs AC=(3, 1) adds one
    count at k=3 and one at k=1.
    """
    hist: dict = {}
    for site in sites:
        for ac in site.get("acs", []):
            hist[ac] = hist.get(ac, 0) + 1
    return hist


def singleton_fraction(hist: dict) -> float:
    """Fraction of alt-allele observations that are singletons (AC=1)."""
    total = sum(hist.values())
    if total == 0:
        return 0.0
    return hist.get(1, 0) / total


def write_sfs_tsv(path: Path, hist: dict) -> None:
    """Persist the SFS histogram as a two-column TSV: ac, n_sites."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("ac\tn_sites\n")
        for k in sorted(hist):
            fh.write(f"{k}\t{hist[k]}\n")
