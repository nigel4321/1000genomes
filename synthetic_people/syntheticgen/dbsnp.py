"""dbSNP rsID injection (M7).

Goal: cohort VCFs need realistic rsID density in the ID column rather
than a sea of ".". The full dbSNP VCF is multi-GB; rather than force
that download, we lean on the fact that the cached ClinVar VCF already
carries dbSNP rs numbers in `INFO/RS`. That's "free" — no extra
download — and gives a few-thousand-record pool per chromosome which is
plenty for injection-density purposes.

Users with a real dbSNP VCF on disk can pass it via `--dbsnp-vcf`; this
module exposes `load_rsid_pool` which accepts either a ClinVar-shaped
VCF (rsIDs in INFO/RS) or a dbSNP-shaped VCF (rsIDs in the ID column).

Injection works the same way as ClinVar injection: pick `density` ×
len(sites) cohort sites at random, overwrite their (pos, ref, alt, id)
with a chosen pool entry, keep the GT block (LD signal) intact, and
re-sort. The simulated REF/ALT bases are random anyway, so swapping
them for plausible dbSNP REF/ALTs is realism-positive.
"""

from __future__ import annotations

import random
import subprocess
from pathlib import Path


# Default fraction of cohort sites to overwrite with a dbSNP-known rsID.
# 0.20 = 20% of records carry an rsID, comfortably above the GA4GH-test
# threshold of "non-empty" while leaving most sites as de-novo coalescent
# variants (no rsID). Tunable via CLI.
DEFAULT_RSID_DENSITY = 0.20


def load_rsid_pool(vcf_path: Path,
                   chromosomes: list[str],
                   max_per_chrom: int | None = 5000,
                   ) -> list[dict]:
    """Stream rsID-bearing variants from a VCF restricted to `chromosomes`.

    Each entry: chrom, pos (int), ref, alt (single string), rsid
    ("rs<digits>"). Multi-allelic sites are split into per-alt rows.
    Symbolic ALTs and length-capped indels are dropped so injected
    records stay writeable through the standard biallelic record path.

    Accepts both ClinVar-style files (rsID in INFO/RS, possibly bare
    digits) and dbSNP-style files (rsID in the ID column with the "rs"
    prefix). Whichever field is non-empty wins; if both are present, the
    ID column is preferred.
    """
    out: list[dict] = []
    for chrom in chromosomes:
        cmd = [
            "bcftools", "query",
            "-r", chrom,
            "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/RS\n",
            str(vcf_path),
        ]
        n_kept = 0
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                if max_per_chrom and n_kept >= max_per_chrom:
                    proc.terminate()
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                c, pos, vid, ref, alt, rs = parts[:6]
                rsid = _normalise_rsid(vid, rs)
                if not rsid:
                    continue
                if len(ref) > 50:
                    continue
                for a in alt.split(","):
                    if a.startswith("<") or len(a) > 50:
                        continue
                    out.append({
                        "chrom": c,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": a,
                        "rsid": rsid,
                    })
                n_kept += 1
    return out


def _normalise_rsid(id_field: str, info_rs: str) -> str:
    """Return a normalised "rsNNN" string or "" if no rsID is present.

    `id_field` may already be `rsNNN` (dbSNP-style) or empty/`.`.
    `info_rs` from ClinVar is bare digits or `.` — wrap with `rs` prefix.
    Multiple rsIDs (semicolon-separated in some files) collapse to the
    first.
    """
    if id_field and id_field != ".":
        first = id_field.split(";")[0].strip()
        if first.startswith("rs"):
            return first
        if first.isdigit():
            return f"rs{first}"
    if info_rs and info_rs != ".":
        first = info_rs.split(",")[0].strip()
        if first.startswith("rs"):
            return first
        if first.isdigit():
            return f"rs{first}"
    return ""


def inject_rsids(sites: list[dict],
                 pool: list[dict],
                 density: float,
                 rng: random.Random,
                 reserve_indices: set[int] | None = None) -> int:
    """Replace `density` × len(sites) cohort sites with dbSNP records.

    The site's GT block (LD signal) is preserved; only coordinates,
    REF/ALT and the ID column are overwritten. `reserve_indices` lets a
    caller exclude sites already claimed by a previous overlay (e.g.
    ClinVar injection) so two overlays don't fight for the same row.

    Sites are sorted by (chrom, pos) on exit. Returns the number of
    injections performed.
    """
    if density <= 0 or not sites or not pool:
        return 0
    chrom_set = {s["chrom"] for s in sites}
    pool_local = [r for r in pool if r["chrom"] in chrom_set]
    if not pool_local:
        return 0

    reserve_indices = reserve_indices or set()
    candidate = [i for i in range(len(sites)) if i not in reserve_indices]
    if not candidate:
        return 0

    n_target = max(1, int(round(density * len(sites))))
    n_target = min(n_target, len(candidate), len(pool_local))

    rng.shuffle(candidate)
    pool_choices = rng.sample(pool_local, n_target)

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
        site["id"] = rec["rsid"]
        injected += 1

    sites.sort(key=lambda s: (s["chrom"], s["pos"]))
    return injected
