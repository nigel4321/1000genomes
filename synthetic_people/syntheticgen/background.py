"""1000G background variant sampling (independent HWE, multi-allelic aware).

Reservoir-samples common variants from local 1000 Genomes VCFs, keeping
per-allele AFs so multi-allelic sites can emit `0|2`, `1|2`, etc. on
draw. Stays as the only background path until the M5 coalescent backbone
lands (at which point this module is gated behind `--legacy-background`).
"""

from __future__ import annotations

import glob
import os
import random
import subprocess
import sys


def load_background_pool(globs: list[str], af_min: float,
                         per_source_limit: int,
                         rng: random.Random) -> list[dict]:
    """Reservoir-sample common variants (AF >= af_min) from each source VCF."""
    sources: list[str] = []
    for g in globs:
        sources.extend(sorted(glob.glob(g)))
    # Deduplicate in case of overlapping globs.
    seen = set()
    sources = [s for s in sources if not (s in seen or seen.add(s))]
    if not sources:
        return []

    pool: list[dict] = []
    for src in sources:
        print(f"  sampling background from {os.path.basename(src)}",
              file=sys.stderr)
        # Pull AF with `-i MAX(INFO/AF)>=af_min` so multi-allelic sites
        # qualify whenever any alt is common enough.
        cmd = [
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
            "-i", f"MAX(INFO/AF)>={af_min}",
            src,
        ]
        reservoir: list[dict] = []
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None
            i = 0
            for line in proc.stdout:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, ref, alt, af_str = parts[:5]
                # Reject symbolic ALTs and length-capped records. Per-allele
                # length cap so a multi-allelic with one long alt doesn't
                # blow up record size.
                alts = alt.split(",")
                if any(a.startswith("<") or len(a) > 50 for a in alts):
                    continue
                if len(ref) > 50:
                    continue
                try:
                    afs = [float(x) for x in af_str.split(",")]
                except ValueError:
                    continue
                if len(afs) != len(alts):
                    continue
                entry = {"chrom": chrom, "pos": int(pos), "id": ".",
                         "ref": ref, "alts": alts, "afs": afs}
                if i < per_source_limit:
                    reservoir.append(entry)
                else:
                    j = rng.randint(0, i)
                    if j < per_source_limit:
                        reservoir[j] = entry
                i += 1
        pool.extend(reservoir)
    return pool


def phased_gt_from_af(af: float, rng: random.Random) -> str:
    """Phased biallelic diploid genotype under Hardy-Weinberg."""
    a = "1" if rng.random() < af else "0"
    b = "1" if rng.random() < af else "0"
    return f"{a}|{b}"


def phased_gt_from_afs(afs: list, rng: random.Random) -> str:
    """Phased diploid genotype for a multi-allelic site.

    `afs` is the list of alt-allele frequencies at the site. Each
    haplotype draws from the categorical distribution {REF, alt_1, ...}
    with P(REF) = 1 - sum(afs). Returns `"a|b"` where `a` and `b` are
    integer allele indices (0 = REF).
    """
    total_alt = sum(afs)
    if total_alt <= 0:
        return "0|0"
    # Cumulative thresholds across alleles 1..k (REF fills the remainder).
    cum = []
    running = 0.0
    for af in afs:
        running += af
        cum.append(running)

    def _draw_allele() -> int:
        r = rng.random()
        for idx, c in enumerate(cum, start=1):
            if r < c:
                return idx
        return 0  # REF

    a = _draw_allele()
    b = _draw_allele()
    return f"{a}|{b}"


def alt_dosages(gt: str, n_alts: int) -> list:
    """Per-alt dosage list: entry k is count of allele (k+1) on the two haplotypes."""
    dosages = [0] * n_alts
    for token in gt.split("|"):
        try:
            idx = int(token)
        except ValueError:
            continue
        if 1 <= idx <= n_alts:
            dosages[idx - 1] += 1
    return dosages


def alt_dosage(gt: str) -> int:
    """Total alt dosage across all alleles — for legacy biallelic call sites."""
    total = 0
    for token in gt.split("|"):
        try:
            idx = int(token)
        except ValueError:
            continue
        if idx >= 1:
            total += 1
    return total


def random_sample_id(rng: random.Random) -> str:
    """HG/NA-prefixed 5-digit ID, mirroring 1000G naming conventions."""
    prefix = rng.choice(("HG", "NA"))
    return f"{prefix}{rng.randint(10000, 99999)}"


def draw_person(candidates: list[dict], background_pool: list[dict],
                n_background: int, rng: random.Random) -> dict:
    """Draw one person's variant set: 1 highlighted + N background.

    Background variants whose drawn genotype turns out hom-ref (0|0) are
    dropped, so the final record count is <= n_background. This mirrors
    real per-sample VCFs which only emit non-reference sites.
    """
    hi = dict(rng.choice(candidates))
    # Highlighted genotype: 70% het, 30% hom-alt, weighted toward the more
    # common "carrier" scenario.
    hi_gt = rng.choices(("0|1", "1|1"), weights=(0.7, 0.3))[0]

    bg_records: list[dict] = []
    if background_pool:
        sample = rng.sample(
            background_pool,
            min(n_background, len(background_pool)),
        )
        for bg in sample:
            gt = phased_gt_from_afs(bg["afs"], rng)
            if alt_dosage(gt) == 0:
                continue
            bg_records.append({**bg, "gt": gt})

    return {
        "sample_id": random_sample_id(rng),
        "highlighted": {**hi, "gt": hi_gt},
        "background": bg_records,
    }
