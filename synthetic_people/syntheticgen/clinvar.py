"""ClinVar fetch + candidate loading + cohort overlay (M7).

Two roles:

* The M1 path: load a small set of pathogenic candidates and use one as
  the per-person "highlighted" variant.
* The M7 path: build a (chrom, pos, ref, alt)-keyed index of ClinVar
  records and either annotate coalescent-produced sites that happen to
  collide with one, or inject a fraction of ClinVar records into the
  cohort so CLNSIG/CLNDN appear at realistic positions. Coalescent sims
  cover positions 1..sim_length while ClinVar sits at real chromosome
  coordinates, so collision-only annotation almost never fires —
  injection is the practical mechanism that makes ClinVar visible in
  the output.

ClinVar's INFO/RS field carries dbSNP rs numbers, so the same cached
file doubles as a rsID source — see `dbsnp.py`.
"""

from __future__ import annotations

import random
import subprocess
import sys
import urllib.request
from pathlib import Path

from .builds import BUILDS


DEFAULT_SIG_FILTER = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
}

# Default fraction of cohort sites to overwrite with injected ClinVar
# records. Roughly the per-genome ClinVar-known-pathogenic density is
# very low; this knob is exposed mainly so the validation suite can see
# CLNSIG-bearing records in every batch.
DEFAULT_CLINVAR_INJECT_DENSITY = 0.01


def _sanitize_info_value(v: str) -> str:
    """Replace characters that would break VCF INFO parsing.

    ClinVar already uses underscores for spaces; this is defensive only.
    """
    return v.replace(";", ",").replace("=", "_").replace(" ", "_")


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
            # Highlighted variants stay single-alt: the "one clinically-
            # highlighted variant per person" concept is inherently one alt.
            if "," in alt or alt.startswith("<") or \
                    len(ref) > 50 or len(alt) > 50:
                continue
            out.append({
                "chrom": chrom,
                "pos": int(pos),
                "id": vid if vid else ".",
                "ref": ref,
                "alts": [alt],
                "afs": [None],  # ClinVar doesn't supply AF; filled by writer.
                "clnsig": _sanitize_info_value(clnsig),
                "clndn": _sanitize_info_value(clndn),
            })
    return out


def load_clinvar_index(clinvar_vcf: Path,
                       chromosomes: list[str],
                       sig_filter: set[str] | None = None,
                       max_per_chrom: int | None = None
                       ) -> list[dict]:
    """Stream ClinVar restricted to `chromosomes` into a flat record list.

    Each entry has chrom/pos/ref/alt/clnsig/clndn/id (the ClinVar VCV id)
    and rsid (from INFO/RS, "" if absent). Multi-allelic sites are kept
    by splitting on commas at parse time. Long indels (>50 bp) and
    symbolic ALTs are dropped — same filter as the highlighted-candidate
    loader, so injected records remain "writeable" through the standard
    record path.

    `sig_filter` defaults to the pathogenic-set the highlighted path
    uses; pass an empty set to skip the CLNSIG filter and load all
    annotated records.
    """
    if sig_filter is None:
        sig_filter = DEFAULT_SIG_FILTER
    out: list[dict] = []
    for chrom in chromosomes:
        cmd = [
            "bcftools", "query",
            "-r", chrom,
            "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t"
                  "%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/RS\n",
            str(clinvar_vcf),
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
                if len(parts) < 8:
                    continue
                c, pos, vid, ref, alt, clnsig, clndn, rs = parts[:8]
                if not clnsig or clnsig == ".":
                    continue
                if sig_filter:
                    sigs = {s.strip() for s in
                            clnsig.replace("|", ",").split(",")}
                    if not sigs & sig_filter:
                        continue
                if len(ref) > 50:
                    continue
                # Multi-allelic sites: emit one row per alt so each is
                # individually injectable as a biallelic record.
                for a in alt.split(","):
                    if a.startswith("<") or len(a) > 50:
                        continue
                    out.append({
                        "chrom": c,
                        "pos": int(pos),
                        "id": vid if vid else ".",
                        "ref": ref,
                        "alt": a,
                        "clnsig": _sanitize_info_value(clnsig),
                        "clndn": _sanitize_info_value(clndn),
                        "rsid": rs if rs and rs != "." else "",
                    })
                n_kept += 1
    return out


def annotate_clinvar(sites: list[dict],
                     clinvar_records: list[dict]) -> int:
    """Overlay CLNSIG/CLNDN/id onto cohort sites that match a ClinVar entry.

    Matches on (chrom, pos, ref, alt[0]); biallelic-only by design (the
    cohort path is biallelic from M5 onwards). Mutates `sites` in place
    and returns the count of annotated sites.
    """
    index: dict = {}
    for r in clinvar_records:
        index[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r
    n = 0
    for s in sites:
        key = (s["chrom"], s["pos"], s["ref"], s["alts"][0])
        rec = index.get(key)
        if rec is None:
            continue
        s["clnsig"] = rec["clnsig"]
        s["clndn"] = rec["clndn"]
        if s.get("id") in (None, "", ".") and rec.get("id"):
            s["id"] = rec["id"]
        n += 1
    return n


def inject_clinvar(sites: list[dict],
                   clinvar_records: list[dict],
                   density: float,
                   rng: random.Random) -> int:
    """Replace `density` × len(sites) cohort sites with ClinVar records.

    The site's GT block (the carrier of LD structure) is preserved; only
    coordinates, REF/ALT, ID and CLNSIG/CLNDN are overwritten. Records
    are picked without replacement from `clinvar_records`, restricted to
    chromosomes that actually appear in `sites`. After injection sites
    remain biallelic and writeable through the standard record path.
    Returns the number of injections performed.

    Sites are sorted by (chrom, pos) on exit, so the cohort site list
    stays monotone for downstream consumers.
    """
    if density <= 0 or not sites or not clinvar_records:
        return 0
    chrom_set = {s["chrom"] for s in sites}
    pool = [r for r in clinvar_records if r["chrom"] in chrom_set]
    if not pool:
        return 0
    n_target = max(1, int(round(density * len(sites))))
    n_target = min(n_target, len(sites), len(pool))

    site_indices = list(range(len(sites)))
    rng.shuffle(site_indices)
    pool_choices = rng.sample(pool, n_target)

    used_keys: set = {(s["chrom"], s["pos"]) for s in sites}
    injected = 0
    cursor = 0
    for rec in pool_choices:
        key = (rec["chrom"], rec["pos"])
        # Skip ClinVar records whose coordinate is already occupied by a
        # cohort site — keeps positions unique without re-deriving them.
        if key in used_keys:
            continue
        # Find the next site index on the same chromosome to swap into.
        target_i = None
        while cursor < len(site_indices):
            i = site_indices[cursor]
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
        site["id"] = rec["id"]
        site["clnsig"] = rec["clnsig"]
        site["clndn"] = rec["clndn"]
        injected += 1

    sites.sort(key=lambda s: (s["chrom"], s["pos"]))
    return injected
