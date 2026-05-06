"""Loader-level tests for the rsID / ClinVar overlays.

The loader functions (`load_rsid_pool`, `load_clinvar_index`,
`load_highlighted_candidates`) shape upstream VCF data into the per-record dicts
the inject_* functions consume. They previously kept records with
ALT="." (or empty ALT) — dbSNP and ClinVar both contain reference-only
entries, especially after batch updates. When such a record was
injected the per-person GT block was preserved, leaving "1|0" calls
against ALT=".", which is invalid VCF: the GT references allele 1 but
allele 1 doesn't exist. Downstream tooling (notably ``bcftools stats``)
crashed on that combination with "Requested allele outside valid
range".

These tests build a small VCF carrying:

  - one ordinary biallelic SNV with an rsID
  - one record where ALT="." with the same rsID

and confirm both loaders drop the malformed record while keeping the
good one.

The unit-level inject_* tests in test_overlays.py work with hand-built
in-memory pools and never see the malformed shape; the bug only
surfaces through the bcftools-driven loader path, which is what this
file exercises.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.clinvar import load_highlighted_candidates, load_clinvar_index
from syntheticgen.dbsnp import load_rsid_pool


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_TABIX = shutil.which("tabix") is not None


_VCF_BODY = """\
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=22,length=50818468,assembly=GRCh38>
##INFO=<ID=RS,Number=.,Type=Integer,Description="dbSNP rs number">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease name">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
22\t1000\trs1000\tA\tG\t.\tPASS\tRS=1000;CLNSIG=Pathogenic;CLNDN=foo
22\t2000\trs2000\tG\t.\t.\tPASS\tRS=2000;CLNSIG=Pathogenic;CLNDN=bar
22\t3000\trs3000\tC\tT\t.\tPASS\tRS=3000;CLNSIG=Pathogenic;CLNDN=baz
"""


@unittest.skipUnless(_HAVE_BCFTOOLS and _HAVE_BGZIP and _HAVE_TABIX,
                     "bcftools / bgzip / tabix not on PATH")
class OverlayLoaderEmptyAltTest(unittest.TestCase):
    """Loaders must drop records with ALT="." so they never reach
    the inject_* path with a payload that produces invalid VCF."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="overlay_loader_")
        plain = Path(cls.tmpdir) / "overlays.vcf"
        plain.write_text(_VCF_BODY)
        subprocess.run(["bgzip", str(plain)], check=True)
        cls.vcf = Path(cls.tmpdir) / "overlays.vcf.gz"
        subprocess.run(["tabix", "-p", "vcf", str(cls.vcf)], check=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_rsid_pool_skips_missing_alt(self):
        pool = load_rsid_pool(self.vcf, ["22"])
        positions = sorted(r["pos"] for r in pool)
        self.assertEqual(
            positions, [1000, 3000],
            msg="loader kept the ALT='.' record at pos 2000; injecting "
                "it would produce 'GT=1|0 against ALT=.' in the cohort "
                "VCF and crash bcftools stats downstream",
        )
        for r in pool:
            self.assertNotIn(r["alt"], ("", "."))

    def test_clinvar_candidates_skips_missing_alt(self):
        # load_highlighted_candidates streams the full file; the fixture
        # VCF is 3 chr22 records, all with CLNSIG=Pathogenic.
        cands = load_highlighted_candidates(self.vcf, {"Pathogenic"})
        positions = sorted(c["pos"] for c in cands)
        self.assertEqual(positions, [1000, 3000])
        for c in cands:
            self.assertTrue(c["alts"][0])
            self.assertNotEqual(c["alts"][0], ".")

    def test_clinvar_index_skips_missing_alt(self):
        idx = load_clinvar_index(self.vcf, ["22"])
        positions = sorted(r["pos"] for r in idx)
        self.assertEqual(positions, [1000, 3000])
        for r in idx:
            self.assertNotIn(r["alt"], ("", "."))


if __name__ == "__main__":
    unittest.main()
