"""ClinVar fetch + candidate loading."""

from __future__ import annotations

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
