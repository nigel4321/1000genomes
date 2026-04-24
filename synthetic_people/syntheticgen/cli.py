"""CLI entry point — wires the package together."""

from __future__ import annotations

import argparse
import random
import shutil
import sys
from pathlib import Path

from .background import load_background_pool, random_sample_id
from .builds import BUILDS
from .clinvar import (
    DEFAULT_SIG_FILTER,
    fetch_clinvar,
    load_highlighted_candidates,
)
from .coalescent import (
    DEFAULT_CHR_LENGTH_MB,
    DEFAULT_DEMO_MODEL,
    DEFAULT_MU,
    DEFAULT_POPULATION,
    DEFAULT_REC_RATE,
    simulate_cohort,
)
from .cohort import draw_cohort_background, person_records_from_cohort
from .sfs import (
    DEFAULT_SFS_ALPHA,
    sfs_histogram,
    singleton_fraction,
    write_sfs_tsv,
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
            "Generate a cohort of synthetic person VCFs. The default path "
            "simulates an LD-correct coalescent with stdpopsim demography; "
            "pass --legacy-background for the M4 1000G-pool + power-law "
            "SFS sampler. Each person receives a ClinVar-highlighted "
            "variant on top of the shared cohort background."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--n", type=int, default=10,
                   help="Cohort size: number of person VCFs to generate")
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
                   help="[legacy] Glob(s) for common-variant source VCFs. "
                        "Pass multiple times to combine sources. "
                        f"Default: {default_bg}")
    p.add_argument("--n-background", type=int, default=500,
                   help="[legacy] Number of shared background sites the "
                        "cohort is drawn over (hom-ref calls dropped "
                        "per person)")
    p.add_argument("--af-min", type=float, default=0.05,
                   help="[legacy] Minimum AF filter when loading the "
                        "coordinate pool from 1000G sources")
    p.add_argument("--sfs-alpha", type=float, default=DEFAULT_SFS_ALPHA,
                   help="[legacy path] Power-law exponent for the cohort "
                        "SFS: P(k) ∝ 1/k^α. α=1.0 is Watterson-neutral; "
                        "α=2.0 (default) biases toward singletons, "
                        "matching gnomAD-like spectra.")
    p.add_argument("--legacy-background", action="store_true",
                   help="Use the M4 1000G-pool + power-law SFS sampler "
                        "instead of the coalescent. No LD, no realistic "
                        "demography — retained for comparison and "
                        "offline-only use.")
    p.add_argument("--chromosomes", default="22",
                   help="[coalescent] Comma-separated chromosomes to "
                        "simulate (e.g. '19,20,21,22'). Shorter list = "
                        "faster run.")
    p.add_argument("--chr-length-mb", type=float,
                   default=DEFAULT_CHR_LENGTH_MB,
                   help="[coalescent] Simulated prefix length per "
                        "chromosome in Mb. 0 = full length (slow on big "
                        "chromosomes, fine for chr22).")
    p.add_argument("--demo-model", default=DEFAULT_DEMO_MODEL,
                   help="[coalescent] stdpopsim demographic model id. "
                        "Pass 'none' for a constant-size Ne=10k "
                        "single-pop msprime draw (no real demography).")
    p.add_argument("--population", default=DEFAULT_POPULATION,
                   help="[coalescent] Sampling population within the demo "
                        "model (e.g. CEU / YRI / CHB for OutOfAfrica_3G09).")
    p.add_argument("--rec-rate", type=float, default=DEFAULT_REC_RATE,
                   help="[coalescent, --demo-model=none only] Uniform "
                        "recombination rate per bp per generation.")
    p.add_argument("--mu", type=float, default=DEFAULT_MU,
                   help="[coalescent, --demo-model=none only] Mutation "
                        "rate per bp per generation.")
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

    if args.legacy_background:
        bg_globs = args.background_glob or [args._default_bg]
        print(f"[legacy] Sampling background coordinates from: {bg_globs}",
              file=sys.stderr)
        background_pool = load_background_pool(
            bg_globs, args.af_min, per_source_limit=5000, rng=rng,
        )
        print(f"  coordinate pool: {len(background_pool)} variants",
              file=sys.stderr)
        if not background_pool:
            print("  [warn] no background sources matched — output VCFs "
                  "will contain only the highlighted variant",
                  file=sys.stderr)
        print(
            f"[legacy] Drawing cohort background: {args.n_background} "
            f"shared sites across {args.n} people (α={args.sfs_alpha})",
            file=sys.stderr,
        )
        cohort_sites = draw_cohort_background(
            background_pool, args.n, args.n_background, args.sfs_alpha,
            rng,
        )
    else:
        chromosomes = [c.strip() for c in args.chromosomes.split(",")
                       if c.strip()]
        demo_model = None if args.demo_model.lower() == "none" \
            else args.demo_model
        print(
            f"Simulating coalescent cohort: {args.n} people across "
            f"chroms {chromosomes} "
            f"(model={demo_model or 'uniform-constant-Ne'}, "
            f"pop={args.population}, length={args.chr_length_mb} Mb)",
            file=sys.stderr,
        )
        cohort_sites = simulate_cohort(
            chromosomes=chromosomes, build=args.build,
            n_people=args.n, length_mb=args.chr_length_mb,
            demo_model=demo_model, population=args.population,
            rec_rate=args.rec_rate, mu=args.mu,
            rng=rng, verbose=True,
        )

    hist = sfs_histogram(cohort_sites)
    total_alts = sum(hist.values())
    singletons = hist.get(1, 0)
    sfrac = singleton_fraction(hist)
    print(
        f"  cohort sites: {len(cohort_sites)}; alt observations: "
        f"{total_alts}; singletons: {singletons} ({sfrac:.1%})",
        file=sys.stderr,
    )

    summary_dir = args.output_dir / "summary"
    sfs_path = summary_dir / "sfs.tsv"
    write_sfs_tsv(sfs_path, hist)
    print(f"  SFS histogram written to {sfs_path}", file=sys.stderr)

    sample_ids = [random_sample_id(rng) for _ in range(args.n)]

    print(f"Writing {args.n} person VCFs into {args.output_dir}",
          file=sys.stderr)
    for i, sid in enumerate(sample_ids):
        hi = dict(rng.choice(candidates))
        hi["gt"] = rng.choices(("0|1", "1|1"), weights=(0.7, 0.3))[0]
        background = person_records_from_cohort(cohort_sites, i)
        person = {
            "sample_id": sid,
            "highlighted": hi,
            "background": background,
        }
        out = args.output_dir / f"person_{i+1:04d}.vcf.gz"
        write_person_vcf(out, person, args.build, rng)
        print(
            f"  [{i+1:>4}/{args.n}] {out.name} — {sid} — "
            f"highlighted {hi['id']} at {hi['chrom']}:{hi['pos']} "
            f"{hi['ref']}>{','.join(hi['alts'])} ({hi['gt']}), "
            f"{len(background)} background records",
            file=sys.stderr,
        )

    print("Done.", file=sys.stderr)
    return 0
