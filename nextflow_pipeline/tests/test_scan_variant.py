"""Tests for bin/scan_variant.py — each status + carrier extraction."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    TARGET_VARIANT,
    above_threshold_variants,
    absent_variants,
    af_unknown_variants,
    below_threshold_variants,
    bin_script,
    default_cohort,
    in_range_variants,
    position_empty_variants,
    require_tools,
    run_script,
    standard_filename,
)
from synthetic_vcf import write_vcf


class ScanVariantTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="scan_variant_test_")

    def _scan(self, vcf_path: str, name: str = "test",
              chrom: str | None = None, pos: int | None = None,
              ref: str | None = None, alt: str | None = None,
              min_af: float = 0.05, max_af: float = 1.0) -> dict:
        out_prefix = os.path.join(self.tmpdir, name)
        run_script(
            bin_script("scan_variant.py"),
            "--vcf", vcf_path,
            "--name", name,
            "--out-prefix", out_prefix,
            "--chrom", chrom or TARGET_VARIANT["chrom"],
            "--pos", str(pos or TARGET_VARIANT["pos"]),
            "--ref", ref or TARGET_VARIANT["ref"],
            "--alt", alt or TARGET_VARIANT["alt"],
            "--min-af", str(min_af),
            "--max-af", str(max_af),
        )
        with open(f"{out_prefix}.variant.json") as fh:
            data = json.load(fh)
        data["_carriers_path"] = f"{out_prefix}.carriers.tsv"
        return data

    def _read_carriers(self, path: str) -> tuple[list[str], list[list[str]]]:
        with open(path) as fh:
            lines = [l for l in fh.read().splitlines() if l]
        header = lines[0].split("\t") if lines else []
        rows = [l.split("\t") for l in lines[1:]]
        return header, rows

    # -------------------------------------------------------------------------
    # Status: not_applicable — chromosome absent from the VCF
    # -------------------------------------------------------------------------

    def test_status_not_applicable_when_chrom_missing(self):
        # Build a chr22-only VCF and query chr15 → NA.
        path = standard_filename(self.tmpdir, "22")
        write_vcf(path, "22", in_range_variants(), default_cohort())
        result = self._scan(path, chrom="15")
        self.assertEqual(result["status"], "not_applicable")
        header, rows = self._read_carriers(result["_carriers_path"])
        self.assertTrue(header)
        self.assertEqual(rows, [])  # header-only, no carriers

    # -------------------------------------------------------------------------
    # Status: position_empty — chrom present, no call at position
    # -------------------------------------------------------------------------

    def test_status_position_empty(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", position_empty_variants(), default_cohort())
        result = self._scan(path)
        self.assertEqual(result["status"], "position_empty")
        _, rows = self._read_carriers(result["_carriers_path"])
        self.assertEqual(rows, [])

    # -------------------------------------------------------------------------
    # Status: absent — position has a record, different REF/ALT
    # -------------------------------------------------------------------------

    def test_status_absent_when_alt_allele_differs(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", absent_variants(), default_cohort())
        result = self._scan(path)  # queries A>G, record is A>C
        self.assertEqual(result["status"], "absent")
        self.assertIn("A>C", result.get("note", ""))

    # -------------------------------------------------------------------------
    # Status: present_in_range — main success path + carriers
    # -------------------------------------------------------------------------

    def test_status_present_in_range_and_carrier_counts(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", in_range_variants(), default_cohort())
        result = self._scan(path, min_af=0.05, max_af=1.0)

        self.assertEqual(result["status"], "present_in_range")
        self.assertEqual(result["variant_id"], "rs12913832")
        self.assertEqual(result["ref"], "A")
        self.assertEqual(result["alt"], "G")
        self.assertEqual(result["an"], 40)  # 20 diploid samples

        # EUR fixed at 1.0 so EUR AF must be 1.0; AFR fixed at 0.0.
        self.assertAlmostEqual(result["eur_af"], 1.0, places=4)
        self.assertAlmostEqual(result["afr_af"], 0.0, places=4)

        # AC integrity: n_heterozygotes + 2*n_homozygous_alt == ac
        self.assertEqual(
            result["ac"],
            result["n_heterozygotes"] + 2 * result["n_homozygous_alt"],
        )

        # Carriers file: all 4 EUR samples are homozygous alt at minimum.
        header, rows = self._read_carriers(result["_carriers_path"])
        self.assertEqual(header[0], "file")
        self.assertEqual(len(rows), result["n_carriers"])
        self.assertGreaterEqual(result["n_carriers"], 4)
        eur_samples = {"HG00096", "HG00097", "HG00099", "HG00100"}
        carrier_samples = {r[1] for r in rows}
        self.assertTrue(eur_samples.issubset(carrier_samples))
        # Each EUR sample should be dosage=2 (homozygous alt).
        eur_rows = [r for r in rows if r[1] in eur_samples]
        for r in eur_rows:
            self.assertEqual(r[8], "2", f"unexpected dosage for {r[1]}: {r[8]}")

    # -------------------------------------------------------------------------
    # Status: present_below_threshold
    # -------------------------------------------------------------------------

    def test_status_present_below_threshold(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", below_threshold_variants(), default_cohort())
        # Only one alt allele in cohort → AF = 1/40 = 0.025, below min 0.05.
        result = self._scan(path, min_af=0.05, max_af=1.0)
        self.assertEqual(result["status"], "present_below_threshold")
        self.assertAlmostEqual(result["af"], 0.025, places=4)
        # Exactly one heterozygous carrier.
        self.assertEqual(result["n_carriers"], 1)
        self.assertEqual(result["n_heterozygotes"], 1)
        self.assertEqual(result["n_homozygous_alt"], 0)

    # -------------------------------------------------------------------------
    # Status: present_above_threshold
    # -------------------------------------------------------------------------

    def test_status_present_above_threshold(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", above_threshold_variants(), default_cohort())
        # AF=1.0 in the cohort; max_af=0.5 → above threshold.
        result = self._scan(path, min_af=0.0, max_af=0.5)
        self.assertEqual(result["status"], "present_above_threshold")
        self.assertAlmostEqual(result["af"], 1.0, places=4)
        # Still emits all carriers (AF gate is file-level only).
        self.assertEqual(result["n_carriers"], 20)
        self.assertEqual(result["n_homozygous_alt"], 20)

    # -------------------------------------------------------------------------
    # Status: present_af_unknown
    # -------------------------------------------------------------------------

    def test_status_present_af_unknown_when_no_af_derivable(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", af_unknown_variants(), default_cohort())
        result = self._scan(path)
        # With AF stripped and all genotypes missing, neither INFO nor
        # fill-tags can produce a value.
        self.assertEqual(result["status"], "present_af_unknown")

    # -------------------------------------------------------------------------
    # chr-prefix resolution: query as "chr15" against a VCF that uses "15"
    # -------------------------------------------------------------------------

    def test_chr_prefix_is_resolved(self):
        path = standard_filename(self.tmpdir, "15")
        write_vcf(path, "15", in_range_variants(), default_cohort())
        result = self._scan(path, chrom="chr15")
        # Pipeline should strip the chr prefix and find the variant.
        self.assertEqual(result["status"], "present_in_range")


if __name__ == "__main__":
    unittest.main()
