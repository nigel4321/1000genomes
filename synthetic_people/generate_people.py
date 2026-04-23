#!/usr/bin/env python3
"""Generate synthetic single-person VCFs.

Each output VCF represents one person with:
  1. one clinically-highlighted variant drawn from ClinVar, annotated with
     its CLNSIG (clinical significance) and CLNDN (disease name);
  2. a realistic background of common variants sampled from local 1000
     Genomes VCFs, with genotypes drawn under Hardy-Weinberg equilibrium
     from each site's cohort allele frequency.

ClinVar is downloaded once and cached under ./cache. Re-runs re-use the
cache. With --seed N the output is deterministic; omit --seed for different
people on every run. Each output VCF is bgzipped + tabix indexed and
compliant with VCFv4.2.

Usage:
    python3 synthetic_people/generate_people.py --n 100 --seed 42
"""

from __future__ import annotations

import argparse
import glob
import os
import random
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


# -----------------------------------------------------------------------------
# Reference builds
# -----------------------------------------------------------------------------

GRCH37_CONTIG_LENGTHS = {
    "1":  249250621, "2":  243199373, "3":  198022430, "4":  191154276,
    "5":  180915260, "6":  171115067, "7":  159138663, "8":  146364022,
    "9":  141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16":  90354753,
    "17":  81195210, "18":  78077248, "19":  59128983, "20":  63025520,
    "21":  48129895, "22":  51304566, "X": 155270560, "Y":  59373566,
    "MT": 16569,
}

GRCH38_CONTIG_LENGTHS = {
    "1":  248956422, "2":  242193529, "3":  198295559, "4":  190214555,
    "5":  181538259, "6":  170805979, "7":  159345973, "8":  145138636,
    "9":  138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16":  90338345,
    "17":  83257441, "18":  80373285, "19":  58617616, "20":  64444167,
    "21":  46709983, "22":  50818468, "X": 156040895, "Y":  57227415,
    "MT": 16569,
}

BUILDS = {
    "GRCh37": {
        "clinvar_url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/"
                       "vcf_GRCh37/clinvar.vcf.gz",
        "reference": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
                     "GCA_000001405.14_GRCh37.p13/"
                     "GCA_000001405.14_GRCh37.p13_genomic.fna.gz",
        "contigs": GRCH37_CONTIG_LENGTHS,
    },
    "GRCh38": {
        "clinvar_url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/"
                       "vcf_GRCh38/clinvar.vcf.gz",
        "reference": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
                     "GCA_000001405.15_GRCh38/"
                     "GCA_000001405.15_GRCh38_genomic.fna.gz",
        "contigs": GRCH38_CONTIG_LENGTHS,
    },
}

DEFAULT_SIG_FILTER = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
}


# -----------------------------------------------------------------------------
# Cache + download
# -----------------------------------------------------------------------------

def _download(url: str, dest: Path) -> None:
    """Stream-download `url` to `dest.part` then rename to `dest`.

    Writing to a .part file means a partial download never masquerades as
    a completed one on the next run.
    """
    tmp = dest.with_name(dest.name + ".part")
    print(f"  downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as resp:
        total = int(resp.getheader("Content-Length") or 0)
        downloaded = 0
        last_report = 0
        with open(tmp, "wb") as fh:
            while True:
                chunk = resp.read(1 << 20)  # 1 MiB
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                if total and downloaded - last_report >= 10 * (1 << 20):
                    pct = downloaded * 100 / total
                    print(f"    {downloaded/1e6:7.1f} / {total/1e6:.1f} MB "
                          f"({pct:.0f}%)", file=sys.stderr)
                    last_report = downloaded
    tmp.rename(dest)


def fetch_clinvar(cache_dir: Path, build: str) -> Path:
    """Ensure the ClinVar VCF + tabix index are cached. Returns VCF path."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    vcf_url = BUILDS[build]["clinvar_url"]
    tbi_url = vcf_url + ".tbi"
    vcf_path = cache_dir / f"clinvar_{build}.vcf.gz"
    tbi_path = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    if not vcf_path.exists():
        _download(vcf_url, vcf_path)
    if not tbi_path.exists():
        _download(tbi_url, tbi_path)
    return vcf_path


# -----------------------------------------------------------------------------
# Variant pools
# -----------------------------------------------------------------------------

def _sanitize_info_value(v: str) -> str:
    """Replace characters that would break VCF INFO parsing."""
    # ClinVar already uses underscores for spaces; this is defensive only.
    return v.replace(";", ",").replace("=", "_").replace(" ", "_")


def load_highlighted_candidates(clinvar_vcf: Path,
                                sig_filter: set[str]) -> list[dict]:
    """Stream ClinVar, keeping records whose CLNSIG matches the filter.

    Returns a list of dicts with chrom/pos/id/ref/alt/clnsig/clndn. Skips
    multi-allelic sites and long indels to keep synthetic output simple.
    """
    cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n",
        str(clinvar_vcf),
    ]
    out: list[dict] = []
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                          stderr=subprocess.DEVNULL) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrom, pos, vid, ref, alt, clnsig, clndn = parts[:7]
            if not clnsig or clnsig == ".":
                continue
            # CLNSIG can be pipe- or comma-joined for multi-condition records.
            sigs = {s.strip() for s in clnsig.replace("|", ",").split(",")}
            if not sigs & sig_filter:
                continue
            # Keep the pool clean: SNV-like, no symbolic ALTs, modest length.
            if "," in alt or alt.startswith("<") or \
                    len(ref) > 50 or len(alt) > 50:
                continue
            out.append({
                "chrom": chrom,
                "pos": int(pos),
                "id": vid if vid else ".",
                "ref": ref,
                "alt": alt,
                "clnsig": _sanitize_info_value(clnsig),
                "clndn": _sanitize_info_value(clndn),
            })
    return out


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


# -----------------------------------------------------------------------------
# Person generation
# -----------------------------------------------------------------------------

def _phased_gt_from_af(af: float, rng: random.Random) -> str:
    """Phased diploid genotype under Hardy-Weinberg."""
    a = "1" if rng.random() < af else "0"
    b = "1" if rng.random() < af else "0"
    return f"{a}|{b}"


def _alt_dosage(gt: str) -> int:
    return sum(1 for c in gt if c == "1")


def _random_sample_id(rng: random.Random) -> str:
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
            gt = _phased_gt_from_af(bg["af"], rng)
            if _alt_dosage(gt) == 0:
                continue
            bg_records.append({**bg, "gt": gt})

    return {
        "sample_id": _random_sample_id(rng),
        "highlighted": {**hi, "gt": hi_gt},
        "background": bg_records,
    }


# -----------------------------------------------------------------------------
# VCF writing
# -----------------------------------------------------------------------------

_HEADER_STATIC = [
    "##fileformat=VCFv4.2",
    "##source=synthetic_people/generate_people.py",
    '##INFO=<ID=AC,Number=A,Type=Integer,'
    'Description="Per-sample alternate allele count">',
    '##INFO=<ID=AN,Number=1,Type=Integer,'
    'Description="Per-sample total allele count">',
    '##INFO=<ID=AF,Number=A,Type=Float,'
    'Description="Per-sample alternate allele frequency (AC/AN)">',
    '##INFO=<ID=HIGHLIGHT,Number=0,Type=Flag,'
    'Description="Clinically-highlighted variant for this synthetic person">',
    '##INFO=<ID=CLNSIG,Number=.,Type=String,'
    'Description="Clinical significance from ClinVar">',
    '##INFO=<ID=CLNDN,Number=.,Type=String,'
    'Description="ClinVar disease name">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


def _contig_sort_key(chrom: str, pos: int,
                     contig_order: dict[str, int]) -> tuple[int, int]:
    """Order chromosomes using the reference dict; unknowns sort last."""
    return (contig_order.get(chrom, len(contig_order)), pos)


def write_person_vcf(out_path: Path, person: dict, build: str) -> Path:
    """Write a single-sample bgzipped+indexed VCF."""
    contigs = BUILDS[build]["contigs"]
    reference = BUILDS[build]["reference"]
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

    with open(plain_path, "w") as fh:
        for line in _HEADER_STATIC[:2]:
            fh.write(line + "\n")
        fh.write(f"##reference={reference}\n")
        for c, length in contigs.items():
            fh.write(f"##contig=<ID={c},length={length}>\n")
        for line in _HEADER_STATIC[2:]:
            fh.write(line + "\n")
        fh.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            f"{person['sample_id']}\n"
        )
        for variant, gt, is_hi in records:
            dosage = _alt_dosage(gt)
            info_parts = [f"AC={dosage}", "AN=2", f"AF={dosage/2:.1f}"]
            if is_hi:
                info_parts.append("HIGHLIGHT")
                if variant.get("clnsig") and variant["clnsig"] != ".":
                    info_parts.append(f"CLNSIG={variant['clnsig']}")
                if variant.get("clndn") and variant["clndn"] != ".":
                    info_parts.append(f"CLNDN={variant['clndn']}")
            fh.write("\t".join([
                variant["chrom"], str(variant["pos"]),
                variant.get("id") or ".",
                variant["ref"], variant["alt"],
                "100", "PASS", ";".join(info_parts), "GT", gt,
            ]) + "\n")

    subprocess.run(["bgzip", "-f", str(plain_path)],
                   check=True, capture_output=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", str_out],
                   check=True, capture_output=True)
    return out_path


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    script_dir = Path(__file__).resolve().parent
    # Default background: 1000G chromosome VCFs sitting alongside this dir.
    default_bg = str(script_dir.parent / "ALL.chr*.phase3_*.genotypes.vcf.gz")

    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
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
    p.add_argument("--build", choices=list(BUILDS), default="GRCh37",
                   help="Reference build; must match background VCFs")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for deterministic output. Omit for "
                        "different people each run.")
    p.add_argument("--background-glob", action="append", default=None,
                   help="Glob(s) for common-variant source VCFs. "
                        "Pass multiple times to combine sources.")
    p.add_argument("--n-background", type=int, default=500,
                   help="Max background variants sampled per person "
                        "(hom-ref draws are dropped before writing)")
    p.add_argument("--af-min", type=float, default=0.05,
                   help="Minimum AF for background variants")
    p.add_argument("--clinvar-sig",
                   default=",".join(sorted(DEFAULT_SIG_FILTER)),
                   help="Comma-separated CLNSIG values to include when "
                        "drawing highlighted variants")
    args = p.parse_args()

    for tool in ("bcftools", "tabix", "bgzip"):
        if not shutil.which(tool):
            sys.exit(f"required tool not on PATH: {tool}")

    rng = random.Random(args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reference build: {args.build}", file=sys.stderr)
    print("Fetching ClinVar (cached across runs)...", file=sys.stderr)
    clinvar_vcf = fetch_clinvar(args.cache_dir, args.build)

    sig_filter = {s.strip() for s in args.clinvar_sig.split(",") if s.strip()}
    print(f"Loading highlighted candidates (CLNSIG in {sorted(sig_filter)})...",
          file=sys.stderr)
    candidates = load_highlighted_candidates(clinvar_vcf, sig_filter)
    print(f"  {len(candidates)} candidates matched", file=sys.stderr)
    if not candidates:
        sys.exit("no ClinVar variants matched the CLNSIG filter — widen "
                 "--clinvar-sig")

    bg_globs = args.background_glob or [default_bg]
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


if __name__ == "__main__":
    sys.exit(main())
