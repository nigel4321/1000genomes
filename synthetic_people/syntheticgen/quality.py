"""Simulated sequencing quality metrics: DP, AD, GQ.

The model here is deliberately simple but distributionally realistic:

* **DP** (read depth) ~ Poisson(λ) per site, where λ is a per-sample
  baseline with modest jitter — so different samples in a cohort have
  different mean depths, and within a sample depth varies site-to-site.
* **AD** (allelic depth) splits DP between REF and ALT per the genotype:
    * `0|0`: all reads REF
    * `1|1`: all reads ALT
    * hets: Binomial(DP, p_alt) with `p_alt ≈ 0.5` biased slightly toward
      REF (~0.475) to match the reference-bias seen in real short-read
      aligners.
  `sum(AD) == DP` always holds, which `bcftools stats` cross-checks.
* **GQ** (genotype quality) is a Phred-like score that degrades when
  depth is low (noisy evidence) or when AD looks inconsistent with the
  called genotype. Floored at 0, capped at 99 per VCF convention.

Everything takes an explicit `random.Random` instance so seeded runs stay
deterministic.
"""

from __future__ import annotations

import math
import random


# Defaults — tunable later via CLI if needed.
DEFAULT_DP_MEAN = 30.0
# Per-sample λ jitter: spread across samples so a 500-person cohort shows
# ~25–35 mean-depth variation, not identical λ for everyone.
DEFAULT_DP_SAMPLE_JITTER_SD = 3.0
# Reference bias on heterozygote sites — ~4–6% ref-leaning is typical for
# Illumina + BWA-MEM short reads. This is the expected ALT-read fraction at
# a het (slightly below 0.5 → the REF allele is over-represented).
HET_ALT_FRAC = 0.475


def sample_lambda(base: float, jitter_sd: float,
                  rng: random.Random) -> float:
    """Per-sample mean depth: base depth + Gaussian jitter, clamped low.

    Clamped at 5 so a pathological jitter draw doesn't produce DP=0 runs.
    """
    lam = base + rng.gauss(0, jitter_sd)
    return max(5.0, lam)


def poisson(lam: float, rng: random.Random) -> int:
    """Knuth's algorithm — fine for the modest λ (5–50) we use here.

    Floors at 0. For larger λ a normal approximation would be faster,
    but we stay well inside the regime where Knuth is exact and fast.
    """
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng.random()
    return k - 1


def _binomial(n: int, p: float, rng: random.Random) -> int:
    """Binomial sampling.

    Uses Python 3.12+ `random.binomialvariate` when available, otherwise
    falls back to n coin flips (fine at our DP range).
    """
    binomfn = getattr(rng, "binomialvariate", None)
    if binomfn is not None:
        return binomfn(n, p)
    return sum(1 for _ in range(n) if rng.random() < p)


def ad_from_gt(gt: str, dp: int, rng: random.Random) -> tuple[int, int]:
    """Split DP into (ref_depth, alt_depth) consistent with the genotype."""
    alt_dosage = sum(1 for c in gt if c == "1")
    if alt_dosage == 0:
        return dp, 0
    if alt_dosage == 2:
        return 0, dp
    # Heterozygote: slight reference bias so hets don't come out perfectly
    # 50/50. Matches what an aligner like BWA-MEM shows in practice.
    alt = _binomial(dp, HET_ALT_FRAC, rng)
    return dp - alt, alt


def gq_from_ad(gt: str, ref_d: int, alt_d: int) -> int:
    """Phred-like genotype quality driven by depth and AD/GT consistency.

    The aim isn't to be GATK-accurate — it's to give a GQ that
      * rises with depth,
      * is high when AD strongly supports the called genotype,
      * drops toward 0 when AD contradicts the call.
    """
    dp = ref_d + alt_d
    if dp == 0:
        return 0
    alt_frac = alt_d / dp

    alt_dosage = sum(1 for c in gt if c == "1")
    if alt_dosage == 0:
        # Homozygous reference: support = P(ref reads) → high when alt_frac≈0.
        support = 1.0 - alt_frac
    elif alt_dosage == 2:
        support = alt_frac
    else:
        # Het: support peaks when alt_frac≈0.5, tails off toward 0 or 1.
        support = 1.0 - 2.0 * abs(alt_frac - 0.5)

    # Phred-like: -10 log10(1 - support). A tight epsilon lets well-supported
    # calls reach the 99 ceiling instead of plateauing at ~40.
    gq = -10.0 * math.log10(max(1.0 - support, 1e-10))
    # Depth cap: low DP can't produce a hugely confident call. Scales up
    # quickly so DP≈30 reaches ~88 and DP≈40+ saturates near 99.
    gq = min(gq, 10.0 * math.log10(max(dp, 1)) * 6.0)
    return max(0, min(99, round(gq)))


def draw_site_quality(gt: str, lam: float,
                      rng: random.Random) -> tuple[int, int, int, int]:
    """Return (DP, ref_AD, alt_AD, GQ) for one record.

    Combines the primitives above; this is what the writer calls per row.
    """
    dp = poisson(lam, rng)
    if dp == 0:
        return 0, 0, 0, 0
    ref_d, alt_d = ad_from_gt(gt, dp, rng)
    gq = gq_from_ad(gt, ref_d, alt_d)
    return dp, ref_d, alt_d, gq
