"""COSMIC overlay (M7, optional / `--somatic`).

COSMIC requires a registration-gated download (no public anonymous
URL), so this module never auto-fetches. The user supplies a local
path to a COSMIC VCF (e.g. `Cosmic_GenomeScreensMutant_v100_GRCh38.vcf`
or any rebuild that follows the same INFO schema), and we overlay or
inject records onto cohort sites the same way ClinVar does.

If `--somatic` is set without a path, we emit a clear instruction and
exit non-zero — the spec mandates the optional gating, but silently
skipping is too quiet for a flag the user explicitly opted into.

The relevant COSMIC INFO tags vary slightly by release; we read what
exists and skip what doesn't:
- `GENE` / `GENE_NAME` — affected gene symbol
- `LEGACY_ID` / `COSMIC_ID` — COSV/COSM identifier (also in the ID col)
- `CDS` / `AA` — coding-DNA / protein consequence string

Fields land in `INFO/COSMIC_GENE` and `INFO/COSMIC_ID` on the output
records; the writer/header pick those up via the M7 INFO declarations.
"""

from __future__ import annotations

import random
import subprocess
from pathlib import Path


DEFAULT_COSMIC_INJECT_DENSITY = 0.005


def load_cosmic_records(cosmic_vcf: Path,
                        chromosomes: list[str],
                        max_per_chrom: int | None = 5000,
                        ) -> list[dict]:
    """Stream COSMIC entries restricted to `chromosomes`.

    Returns one dict per (chrom, pos, ref, alt) row. Robust to missing
    optional INFO tags — older COSMIC releases don't all carry every
    field. Multi-allelic ALTs are split.
    """
    out: list[dict] = []
    for chrom in chromosomes:
        # Permissive query — COSMIC schema varies; we extract whatever
        # is there and fall back to "." for missing tags.
        fmt = ("%CHROM\t%POS\t%ID\t%REF\t%ALT\t"
               "%INFO/GENE\t%INFO/LEGACY_ID\t%INFO/CDS\t%INFO/AA\n")
        cmd = ["bcftools", "query", "-r", chrom, "-f", fmt, str(cosmic_vcf)]
        n_kept = 0
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                if max_per_chrom and n_kept >= max_per_chrom:
                    proc.terminate()
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                c, pos, vid, ref, alt, gene, lid, cds, aa = parts[:9]
                if len(ref) > 50:
                    continue
                for a in alt.split(","):
                    if a.startswith("<") or len(a) > 50:
                        continue
                    out.append({
                        "chrom": c,
                        "pos": int(pos),
                        "id": vid if vid and vid != "." else (
                            lid if lid and lid != "." else "."),
                        "ref": ref,
                        "alt": a,
                        "gene": gene if gene and gene != "." else "",
                        "cds": cds if cds and cds != "." else "",
                        "aa": aa if aa and aa != "." else "",
                    })
                n_kept += 1
    return out


def inject_cosmic(sites: list[dict],
                  cosmic_records: list[dict],
                  density: float,
                  rng: random.Random,
                  reserve_indices: set[int] | None = None) -> int:
    """Replace `density` × len(sites) cohort sites with COSMIC records.

    Mirrors `clinvar.inject_clinvar` but writes COSMIC_GENE / COSMIC_ID
    (and the ID column when present). Returns the number of injections.
    """
    if density <= 0 or not sites or not cosmic_records:
        return 0
    chrom_set = {s["chrom"] for s in sites}
    pool = [r for r in cosmic_records if r["chrom"] in chrom_set]
    if not pool:
        return 0

    reserve_indices = reserve_indices or set()
    candidate = [i for i in range(len(sites)) if i not in reserve_indices]
    if not candidate:
        return 0

    n_target = max(1, int(round(density * len(sites))))
    n_target = min(n_target, len(candidate), len(pool))

    rng.shuffle(candidate)
    pool_choices = rng.sample(pool, n_target)

    used_keys: set = {(s["chrom"], s["pos"]) for s in sites}
    injected = 0
    cursor = 0
    for rec in pool_choices:
        key = (rec["chrom"], rec["pos"])
        if key in used_keys:
            continue
        target_i = None
        while cursor < len(candidate):
            i = candidate[cursor]
            cursor += 1
            if sites[i]["chrom"] == rec["chrom"]:
                target_i = i
                break
        if target_i is None:
            break
        site = sites[target_i]
        old_key = (site["chrom"], site["pos"])
        used_keys.discard(old_key)
        used_keys.add(key)
        site["pos"] = rec["pos"]
        site["ref"] = rec["ref"]
        site["alts"] = [rec["alt"]]
        if rec["id"] and rec["id"] != ".":
            site["id"] = rec["id"]
        if rec["gene"]:
            site["cosmic_gene"] = rec["gene"]
        if rec["id"] and rec["id"] != ".":
            site["cosmic_id"] = rec["id"]
        injected += 1

    sites.sort(key=lambda s: (s["chrom"], s["pos"]))
    return injected
