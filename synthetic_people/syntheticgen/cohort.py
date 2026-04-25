"""Cohort-level site generation — shared coordinates, per-person genotypes.

From M4 onwards the generator simulates an N-sample cohort as a single
event: one pass selects the variable sites and their allele frequencies,
and each site is then "populated" by assigning alt alleles into specific
haplotype slots. Diploid GTs per person fall out of pairing consecutive
slot assignments. This replaces the M1–M3 per-person independent sampler.

Exact slot assignment (sampling without replacement from 2n haplotypes)
preserves the drawn minor allele count exactly, so the realised cohort
SFS matches the SFS we drew from — no smoothing from independent HWE
resampling at every site. The pairing step inherits a mild
hypergeometric-vs-HWE correction, which is if anything more realistic
for a finite cohort than strict HWE draws.
"""

from __future__ import annotations

import random

from .sfs import DEFAULT_SFS_ALPHA, draw_allele_counts


def assign_haplotypes(n_haplotypes: int, allele_counts: list,
                      rng: random.Random) -> list:
    """Assign alt alleles to haplotype slots.

    Returns a length-n_haplotypes list where each entry is an allele
    index (0 = REF, 1..k = alt_1..alt_k). The assignment is uniformly
    random across slots without replacement, so the realised count of
    each allele equals the input `allele_counts` exactly.
    """
    total_alt = sum(allele_counts)
    if total_alt > n_haplotypes:
        raise ValueError(
            f"total alt count {total_alt} exceeds haplotype slots "
            f"{n_haplotypes}"
        )
    slots = [0] * n_haplotypes
    positions = list(range(n_haplotypes))
    rng.shuffle(positions)
    cursor = 0
    for alt_index, count in enumerate(allele_counts, start=1):
        for _ in range(count):
            slots[positions[cursor]] = alt_index
            cursor += 1
    return slots


def _gts_from_slots(slots: list, n_people: int) -> list:
    """Pair consecutive haplotype slots into phased diploid GT strings."""
    return [f"{slots[2 * i]}|{slots[2 * i + 1]}" for i in range(n_people)]


def draw_cohort_background(pool: list, n_people: int, n_sites: int,
                           alpha: float,
                           rng: random.Random) -> list:
    """Generate the cohort's shared background sites.

    `pool` provides candidate coordinates (chrom, pos, ref, alts, ...).
    For each chosen site, allele counts are redrawn from the power-law
    SFS (the 1000G source AFs are ignored — only its genomic coordinates
    and allele strings carry through). Alt alleles are then placed into
    specific haplotype slots so that each person's genotype at the site
    is consistent with every other person in the cohort.
    """
    if not pool or n_people < 1 or n_sites < 1:
        return []
    n_haplotypes = 2 * n_people
    sample_n = min(n_sites, len(pool))
    chosen = rng.sample(pool, sample_n)
    sites = []
    for entry in chosen:
        n_alts = len(entry["alts"])
        counts = draw_allele_counts(n_haplotypes, n_alts, alpha, rng)
        afs = [c / n_haplotypes for c in counts]
        slots = assign_haplotypes(n_haplotypes, counts, rng)
        gts = _gts_from_slots(slots, n_people)
        sites.append({
            "chrom": entry["chrom"],
            "pos": entry["pos"],
            "id": entry.get("id", "."),
            "ref": entry["ref"],
            "alts": list(entry["alts"]),
            "afs": afs,
            "acs": counts,
            "gts": gts,
        })
    return sites


def person_records_from_cohort(sites: list, person_index: int) -> list:
    """Project cohort sites down to one person's non-hom-ref records.

    Carries the M7 annotation fields (clnsig, clndn, cosmic_id,
    cosmic_gene) through to the per-person record dict so the writer can
    emit the corresponding INFO tags for that sample.
    """
    records = []
    for site in sites:
        gt = site["gts"][person_index]
        if all(tok == "0" for tok in gt.split("|")):
            continue
        rec = {
            "chrom": site["chrom"],
            "pos": site["pos"],
            "id": site["id"],
            "ref": site["ref"],
            "alts": site["alts"],
            "afs": site["afs"],
            "gt": gt,
        }
        for key in ("clnsig", "clndn", "cosmic_id", "cosmic_gene"):
            if site.get(key):
                rec[key] = site[key]
        records.append(rec)
    return records
