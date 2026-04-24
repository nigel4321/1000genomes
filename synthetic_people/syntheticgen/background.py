"""Legacy 1000G background variant sampling (independent HWE).

This is the M0 behaviour — kept for parity while coalescent-based generation
(M5+) is being developed. CLI exposes it via `--legacy-background` eventually;
for M1 it is still the only path.
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
        cmd = [
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
            "-i", f"INFO/AF>={af_min}",
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
                if "," in alt or alt.startswith("<") or \
                        len(ref) > 50 or len(alt) > 50:
                    continue
                try:
                    af = float(af_str)
                except ValueError:
                    continue
                entry = {"chrom": chrom, "pos": int(pos), "id": ".",
                         "ref": ref, "alt": alt, "af": af}
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
    """Phased diploid genotype under Hardy-Weinberg."""
    a = "1" if rng.random() < af else "0"
    b = "1" if rng.random() < af else "0"
    return f"{a}|{b}"


def alt_dosage(gt: str) -> int:
    return sum(1 for c in gt if c == "1")


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
            gt = phased_gt_from_af(bg["af"], rng)
            if alt_dosage(gt) == 0:
                continue
            bg_records.append({**bg, "gt": gt})

    return {
        "sample_id": random_sample_id(rng),
        "highlighted": {**hi, "gt": hi_gt},
        "background": bg_records,
    }
