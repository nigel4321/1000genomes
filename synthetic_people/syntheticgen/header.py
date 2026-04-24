"""VCF header assembly — contig lines with assembly=, INFO, FORMAT, ALT."""

from __future__ import annotations

from .builds import BUILDS


# Core INFO lines used since M1.
_INFO_CORE = [
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
]

# Structural-variant INFO lines declared now; values populated in M8.
_INFO_SV = [
    '##INFO=<ID=SVTYPE,Number=1,Type=String,'
    'Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,'
    'Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,'
    'Description="End position of the variant described in this record">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,'
    'Description="Confidence interval around POS for imprecise variants">',
]

# Symbolic ALT declarations for SVs (M8 populates these records).
_ALT_SV = [
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=INS,Description="Insertion">',
]

# FORMAT lines. GT has values from M1; DP/GQ/AD become meaningful in M2.
_FORMAT = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,'
    'Description="Read depth at this position (simulated)">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,'
    'Description="Genotype quality (simulated)">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,'
    'Description="Allelic depths for the ref and alt alleles (simulated)">',
]


def build_header(build: str, sample_id: str) -> str:
    """Return the full VCF header for a single-sample file.

    Per M1, the header fully declares the tag set the project will use
    across milestones: `GT`, `DP`, `GQ`, `AD` FORMAT and `SVTYPE`, `SVLEN`,
    `END`, `CIPOS` INFO plus the `<DEL>`/`<DUP>`/`<INV>`/`<INS>` symbolic
    ALTs. Records before M2/M8 only fill in `GT`; downstream tools tolerate
    declared-but-unused tags, and having them declared early means the
    header stops changing shape as later milestones land.
    """
    info = BUILDS[build]
    assembly = info["assembly"]
    reference = info["reference"]
    contigs = info["contigs"]

    lines: list[str] = [
        "##fileformat=VCFv4.2",
        "##source=synthetic_people/generate_people.py",
        f"##reference={reference}",
    ]
    for chrom, length in contigs.items():
        lines.append(
            f"##contig=<ID={chrom},length={length},assembly={assembly}>"
        )
    lines.extend(_INFO_CORE)
    lines.extend(_INFO_SV)
    lines.extend(_ALT_SV)
    lines.extend(_FORMAT)
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id
    )
    return "\n".join(lines) + "\n"
