"""Ti/Tv calibrator — bias SNV alt-allele selection toward a target ratio.

Transitions (Ti) are purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T)
substitutions. Transversions (Tv) are everything else. Unbiased uniform
drawing of an ALT base from {A,C,G,T}\\{REF} produces Ti:Tv = 4:8 = 0.5,
which is nowhere near the empirical human WGS ratio of ~2.1 (Ti/Tv).

This module exists so the coalescent / de-novo SNV path (M5 onwards) can
hit the realistic ratio without any post-hoc rejection sampling. The
current 1000G-backed M1-M4 sampler doesn't need it (real data is
naturally Ti-enriched).
"""

from __future__ import annotations

import random


DEFAULT_TARGET_TITV = 2.1

# Purine↔purine and pyrimidine↔pyrimidine pairs are transitions.
_TRANSITION_PARTNER = {"A": "G", "G": "A", "C": "T", "T": "C"}
_BASES: tuple = ("A", "C", "G", "T")


def is_transition(ref: str, alt: str) -> bool:
    """True if the single-base substitution REF→ALT is a transition."""
    return _TRANSITION_PARTNER.get(ref.upper()) == alt.upper()


def choose_alt(ref: str, rng: random.Random,
               target: float = DEFAULT_TARGET_TITV) -> str | None:
    """Sample a non-REF base so long-run Ti/Tv converges on `target`.

    For a given REF, exactly one of the three candidate ALT bases is a
    transition and the other two are transversions. To hit a target
    Ti/Tv of `r`, we give the transition allele weight `r` and each
    transversion allele weight `0.5` (so total Tv weight is 1). The
    asymptotic Ti-to-Tv acceptance ratio is then `r : 1`.

    Returns None if `ref` is not a standard base (N, or multi-base REF
    from an indel — callers should route indels elsewhere).
    """
    if target <= 0:
        raise ValueError(f"target Ti/Tv must be positive, got {target}")
    ref = ref.upper()
    if ref not in _BASES:
        return None
    ti_partner = _TRANSITION_PARTNER[ref]
    alts = [b for b in _BASES if b != ref]
    weights = [target if b == ti_partner else 0.5 for b in alts]
    return rng.choices(alts, weights=weights, k=1)[0]


def titv_ratio(pairs) -> float:
    """Compute Ti/Tv on an iterable of (ref, alt) pairs.

    Pairs where either allele is not a single standard base are skipped
    — indels, multi-allelic records, and N-containing refs don't count.
    Returns `float('inf')` if the stream has transitions but no
    transversions, and 0.0 if there are no transitions at all.
    """
    ti = 0
    tv = 0
    for ref, alt in pairs:
        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref.upper() not in _BASES or alt.upper() not in _BASES:
            continue
        if ref.upper() == alt.upper():
            continue
        if is_transition(ref, alt):
            ti += 1
        else:
            tv += 1
    if tv == 0:
        return float("inf") if ti else 0.0
    return ti / tv
