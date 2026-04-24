"""CLI entry point — wires the package together."""

from __future__ import annotations

import argparse
import random
import shutil
import sys
from pathlib import Path

from .background import draw_person, load_background_pool
from .builds import BUILDS
from .clinvar import (
    DEFAULT_SIG_FILTER,
    fetch_clinvar,
    load_highlighted_candidates,
)
from .writer import write_person_vcf


OPTIONAL_PYTHON_DEPS = [
    # These aren't used in M1 but will be required from M2+. `--check-deps`
    # reports their absence so a user can install ahead of time.
    ("numpy", "M2+ (quality metrics, sampling)"),
    ("msprime", "M5+ (coalescent simulation)"),
    ("tskit", "M5+ (tree-sequence handling)"),
    ("stdpopsim", "M6+ (human demographic models)"),
    ("matplotlib", "M10 (validation plots)"),
    ("allel", "M10 (scikit-allel PCA / LD decay)"),
]


def _check_deps(verbose: bool = True) -> int:
    """Check required binaries and (optionally) Python deps. Returns exit code."""
    missing_bins: list[str] = []
    for tool in ("bcftools", "tabix", "bgzip"):
        if not shutil.which(tool):
            missing_bins.append(tool)

    missing_py: list[tuple[str, str]] = []
    for mod, reason in OPTIONAL_PYTHON_DEPS:
        try:
            __import__(mod)
        except ImportError:
            missing_py.append((mod, reason))

    if verbose:
        if not missing_bins:
            print("htslib binaries: OK (bcftools, tabix, bgzip)",
                  file=sys.stderr)
        else:
            print("htslib binaries MISSING: " + ", ".join(missing_bins),
                  file=sys.stderr)
        if not missing_py:
            print("optional Python deps: all present", file=sys.stderr)
        else:
            print("optional Python deps (not required for M1):",
                  file=sys.stderr)
            for mod, reason in missing_py:
                print(f"  - {mod:<12} needed for {reason}", file=sys.stderr)
    # Only hard-fail on missing binaries. Python deps are optional in M1.
    return 1 if missing_bins else 0


def _parser(script_dir: Path) -> argparse.ArgumentParser:
    default_bg = str(script_dir.parent / "ALL.chr*.phase3_*.genotypes.vcf.gz")

    p = argparse.ArgumentParser(
        description=(
            "Generate synthetic single-person VCFs with ClinVar-highlighted "
            "variants and an HWE-sampled 1000 Genomes background."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--n", type=int, default=10,
                   help="Number of person VCFs to generate")
    p.add_argument("--output-dir", type=Path,
                   default=script_dir / "out",
                   help="Where to write person_<N>.vcf.gz")
    p.add_argument("--cache-dir", type=Path,
                   default=script_dir / "cache",
                   help="Where ClinVar is downloaded and cached")
    p.add_argument("--build", choices=list(BUILDS), default="GRCh38",
                   help="Reference build; must match background VCFs")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for deterministic output. Omit for "
                        "different people each run.")
    p.add_argument("--background-glob", action="append", default=None,
                   help="Glob(s) for common-variant source VCFs. "
                        "Pass multiple times to combine sources. "
                        f"Default: {default_bg}")
    p.add_argument("--n-background", type=int, default=500,
                   help="Max background variants sampled per person "
                        "(hom-ref draws are dropped before writing)")
    p.add_argument("--af-min", type=float, default=0.05,
                   help="Minimum AF for background variants")
    p.add_argument("--clinvar-sig",
                   default=",".join(sorted(DEFAULT_SIG_FILTER)),
                   help="Comma-separated CLNSIG values to include when "
                        "drawing highlighted variants")
    p.add_argument("--check-deps", action="store_true",
                   help="Check htslib binaries and optional Python deps, "
                        "then exit")
    p.set_defaults(_default_bg=default_bg)
    return p


def main(argv: list[str] | None = None) -> int:
    script_dir = Path(__file__).resolve().parent.parent
    args = _parser(script_dir).parse_args(argv)

    if args.check_deps:
        return _check_deps()

    # Hard-fail on missing htslib binaries even without --check-deps.
    for tool in ("bcftools", "tabix", "bgzip"):
        if not shutil.which(tool):
            sys.exit(f"required tool not on PATH: {tool}")

    rng = random.Random(args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reference build: {args.build}", file=sys.stderr)
    print("Fetching ClinVar (cached across runs)...", file=sys.stderr)
    clinvar_vcf = fetch_clinvar(args.cache_dir, args.build)

    sig_filter = {s.strip() for s in args.clinvar_sig.split(",") if s.strip()}
    print(
        f"Loading highlighted candidates (CLNSIG in {sorted(sig_filter)})...",
        file=sys.stderr,
    )
    candidates = load_highlighted_candidates(clinvar_vcf, sig_filter)
    print(f"  {len(candidates)} candidates matched", file=sys.stderr)
    if not candidates:
        sys.exit("no ClinVar variants matched the CLNSIG filter — widen "
                 "--clinvar-sig")

    bg_globs = args.background_glob or [args._default_bg]
    print(f"Sampling background from: {bg_globs}", file=sys.stderr)
    background_pool = load_background_pool(
        bg_globs, args.af_min, per_source_limit=5000, rng=rng,
    )
    print(f"  background pool: {len(background_pool)} variants",
          file=sys.stderr)
    if not background_pool:
        print("  [warn] no background sources matched — output VCFs will "
              "contain only the highlighted variant", file=sys.stderr)

    print(f"Generating {args.n} person VCFs into {args.output_dir}",
          file=sys.stderr)
    for i in range(args.n):
        person = draw_person(candidates, background_pool,
                             args.n_background, rng)
        out = args.output_dir / f"person_{i+1:04d}.vcf.gz"
        write_person_vcf(out, person, args.build)
        hi = person["highlighted"]
        print(
            f"  [{i+1:>4}/{args.n}] {out.name} — {person['sample_id']} — "
            f"highlighted {hi['id']} at {hi['chrom']}:{hi['pos']} "
            f"{hi['ref']}>{hi['alt']} ({hi['gt']}), "
            f"{len(person['background'])} background records",
            file=sys.stderr,
        )

    print("Done.", file=sys.stderr)
    return 0
