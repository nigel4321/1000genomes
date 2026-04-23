"""Generate synthetic VCF files that mimic 1000 Genomes Phase 3 data.

Output is bgzipped + tabix-indexed and matches the conventions the pipeline
expects (VCFv4.1 header, ##reference, ##source=1000GenomesPhase3Pipeline,
per-super-pop AF INFO fields, phased diploid genotypes, filename stamped with
`phase3_shapeit2_mvncall_integrated_v5b.<date>`).

Stdlib only — no pysam/cyvcf2 dependency.
"""

from __future__ import annotations

import os
import random
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Iterable


# Real 1000G sample IDs grouped by super-population. Using real IDs keeps the
# synthetic data recognisable to anyone who works with this dataset.
SUPERPOPULATIONS = {
    "EUR": ["HG00096", "HG00097", "HG00099", "HG00100"],
    "AFR": ["NA18486", "NA18488", "NA18489", "NA18498"],
    "EAS": ["HG00403", "HG00404", "HG00406", "HG00407"],
    "SAS": ["HG03006", "HG03007", "HG03008", "HG03009"],
    "AMR": ["HG01112", "HG01113", "HG01119", "HG01121"],
}

DEFAULT_SAMPLES: list[str] = [s for pop in SUPERPOPULATIONS.values() for s in pop]

GRCH37_CONTIG_LENGTHS = {
    "1":  249250621, "2":  243199373, "3":  198022430, "4":  191154276,
    "5":  180915260, "6":  171115067, "7":  159138663, "8":  146364022,
    "9":  141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16":  90354753,
    "17":  81195210, "18":  78077248, "19":  59128983, "20":  63025520,
    "21":  48129895, "22":  51304566, "X": 155270560, "Y":  59373566,
}


@dataclass
class Variant:
    """A single site to emit in the synthetic VCF.

    Either provide `genotypes` directly (one string per sample, e.g. "0|1")
    or provide `af_by_pop` and let the generator draw phased genotypes.
    """

    pos: int
    ref: str
    alt: str
    variant_id: str = "."
    # Direct-specification path: one genotype string per sample.
    genotypes: list[str] | None = None
    # Sampling path: per-super-pop alt-allele frequency. Unseeded callers
    # should construct the generator with a fixed seed for reproducibility.
    af_by_pop: dict[str, float] | None = None
    # If True, omit the AF INFO field entirely (to exercise recomputation /
    # `present_af_unknown` paths in the pipeline). AC/AN are still emitted.
    strip_af_info: bool = False
    # Force all genotypes to "./." — used to exercise AF-unknown edge cases.
    all_missing: bool = False
    extra_info: dict[str, str] = field(default_factory=dict)


@dataclass
class SyntheticCohort:
    samples: list[str] = field(default_factory=lambda: list(DEFAULT_SAMPLES))
    # Which super-pop each sample belongs to. Samples not listed here default
    # to "EUR" for computation convenience.
    pop_by_sample: dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        if not self.pop_by_sample:
            self.pop_by_sample = {
                s: pop
                for pop, members in SUPERPOPULATIONS.items()
                for s in members
            }

    def samples_in_pop(self, pop: str) -> list[str]:
        return [s for s in self.samples if self.pop_by_sample.get(s) == pop]


def _draw_phased_gt(rng: random.Random, af: float) -> str:
    """Draw a phased diploid genotype given an alt-allele frequency."""
    a = "1" if rng.random() < af else "0"
    b = "1" if rng.random() < af else "0"
    return f"{a}|{b}"


def _realize_genotypes(variant: Variant, cohort: SyntheticCohort,
                       rng: random.Random) -> list[str]:
    if variant.all_missing:
        return ["./." for _ in cohort.samples]
    if variant.genotypes is not None:
        if len(variant.genotypes) != len(cohort.samples):
            raise ValueError(
                f"variant at pos {variant.pos}: got {len(variant.genotypes)} "
                f"genotypes for {len(cohort.samples)} samples"
            )
        return list(variant.genotypes)
    af_by_pop = variant.af_by_pop or {}
    default_af = af_by_pop.get("ALL", 0.0)
    out: list[str] = []
    for s in cohort.samples:
        pop = cohort.pop_by_sample.get(s, "EUR")
        af = af_by_pop.get(pop, default_af)
        out.append(_draw_phased_gt(rng, af))
    return out


def _compute_pop_stats(gts: list[str], samples: list[str],
                       cohort: SyntheticCohort) -> tuple[int, int, dict[str, float]]:
    """Return (AC, AN, per-pop AF)."""
    ac = an = 0
    per_pop_ac: dict[str, int] = {}
    per_pop_an: dict[str, int] = {}
    for s, gt in zip(samples, gts):
        pop = cohort.pop_by_sample.get(s, "EUR")
        for allele in gt.replace("|", "/").split("/"):
            if allele in ("0", "1"):
                an += 1
                per_pop_an[pop] = per_pop_an.get(pop, 0) + 1
                if allele == "1":
                    ac += 1
                    per_pop_ac[pop] = per_pop_ac.get(pop, 0) + 1
    pop_af = {}
    for pop in ("EAS", "AMR", "AFR", "EUR", "SAS"):
        pac = per_pop_ac.get(pop, 0)
        pan = per_pop_an.get(pop, 0)
        pop_af[pop] = (pac / pan) if pan else 0.0
    return ac, an, pop_af


def _format_info(ac: int, an: int, pop_af: dict[str, float],
                 variant: Variant) -> str:
    af_overall = (ac / an) if an else 0.0
    pieces: list[str] = []
    if not variant.strip_af_info:
        pieces.append(f"AC={ac}")
        pieces.append(f"AF={af_overall:.6g}")
    pieces.append(f"AN={an}")
    pieces.append(f"NS={an // 2 if an else 0}")
    if not variant.strip_af_info:
        # Order deliberately matches real 1000G INFO line ordering.
        for pop in ("EAS", "AMR", "AFR", "EUR", "SAS"):
            pieces.append(f"{pop}_AF={pop_af.get(pop, 0.0):.6g}")
    pieces.append("VT=SNP")
    for k, v in variant.extra_info.items():
        pieces.append(f"{k}={v}")
    return ";".join(pieces)


def _build_header(chrom: str, cohort: SyntheticCohort,
                  file_date: str) -> list[str]:
    lines = [
        "##fileformat=VCFv4.1",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        f"##fileDate={file_date}",
        "##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/"
        "reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
        "##source=1000GenomesPhase3Pipeline",
    ]
    # Full GRCh37 contig set keeps the header representative and means the
    # pipeline's chr-resolution logic (bare '15' vs 'chr15') is exercised.
    for c, length in GRCH37_CONTIG_LENGTHS.items():
        lines.append(f"##contig=<ID={c},assembly=b37,length={length}>")
    lines += [
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of '
        'alternate alleles in called genotypes">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele '
        'frequency in the range (0,1)">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of '
        'alleles in called genotypes">',
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples '
        'with data">',
        '##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency '
        'in the EAS populations calculated from AC and AN">',
        '##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency '
        'in the EUR populations calculated from AC and AN">',
        '##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency '
        'in the AFR populations calculated from AC and AN">',
        '##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency '
        'in the AMR populations calculated from AC and AN">',
        '##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency '
        'in the SAS populations calculated from AC and AN">',
        '##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type '
        'of variant the line represents">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(cohort.samples)
    )
    return lines


def write_vcf(
    output_path: str,
    chrom: str,
    variants: Iterable[Variant],
    cohort: SyntheticCohort | None = None,
    file_date: str = "20130502",
    seed: int = 1000,
) -> str:
    """Write a bgzipped, tabix-indexed VCF. Returns the .vcf.gz path.

    `output_path` should end in `.vcf.gz`. A matching `.tbi` is generated
    alongside it.
    """
    if not output_path.endswith(".vcf.gz"):
        raise ValueError("output_path must end in .vcf.gz")
    if not shutil.which("bgzip") or not shutil.which("tabix"):
        raise RuntimeError("bgzip and tabix are required to build test VCFs")

    cohort = cohort or SyntheticCohort()
    rng = random.Random(seed)

    # Variants must be emitted in sorted order for tabix to index correctly.
    vs = sorted(variants, key=lambda v: v.pos)

    plain_path = output_path[:-len(".gz")]  # .vcf
    with open(plain_path, "w") as fh:
        for line in _build_header(chrom, cohort, file_date):
            fh.write(line + "\n")
        for v in vs:
            gts = _realize_genotypes(v, cohort, rng)
            ac, an, pop_af = _compute_pop_stats(gts, cohort.samples, cohort)
            info = _format_info(ac, an, pop_af, v)
            fields = [
                chrom, str(v.pos), v.variant_id, v.ref, v.alt,
                "100", "PASS", info, "GT", *gts,
            ]
            fh.write("\t".join(fields) + "\n")

    # bgzip overwrites .vcf.gz if it exists (-f), leaves no .vcf behind.
    subprocess.run(["bgzip", "-f", plain_path], check=True,
                   capture_output=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", output_path], check=True,
                   capture_output=True)
    return output_path


def standard_filename(dir_path: str, chrom: str,
                      date_stamp: str = "20130502") -> str:
    """Return a pipeline-representative filename inside `dir_path`."""
    base = (
        f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
        f"{date_stamp}.genotypes.vcf.gz"
    )
    return os.path.join(dir_path, base)
