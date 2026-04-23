"""Tests for bin/inspect_vcf.py — metadata collection from a VCF."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    bin_script,
    build_standard_cohort_vcf,
    default_cohort,
    require_tools,
    run_script,
    standard_filename,
)
from synthetic_vcf import SyntheticCohort, Variant, write_vcf


class InspectVcfTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="inspect_vcf_test_")
        cls.vcf_path = build_standard_cohort_vcf(
            standard_filename(cls.tmpdir, "15"), chrom="15"
        )

    def _run(self, vcf_path: str, name: str) -> dict:
        proc = run_script(bin_script("inspect_vcf.py"), vcf_path, name)
        return json.loads(proc.stdout)

    def test_emits_expected_top_level_fields(self):
        meta = self._run(self.vcf_path, "chr15_synth")
        expected = {
            "name", "file", "contigs", "n_variants", "n_samples",
            "sample_hash", "reference", "size_bytes", "pipeline_tag",
            "date_stamp",
        }
        self.assertTrue(expected.issubset(meta.keys()))
        self.assertEqual(meta["name"], "chr15_synth")

    def test_sample_count_matches_cohort(self):
        meta = self._run(self.vcf_path, "chr15_synth")
        self.assertEqual(meta["n_samples"], len(default_cohort().samples))

    def test_contigs_match_emitted_chromosome(self):
        # The indexed contig is the one actually seen in the body, not every
        # ##contig line in the header — tabix only records contigs with data.
        meta = self._run(self.vcf_path, "chr15_synth")
        self.assertEqual(meta["contigs"], ["15"])

    def test_variant_count(self):
        meta = self._run(self.vcf_path, "chr15_synth")
        # Standard cohort emits 3 variants in in_range_variants().
        self.assertEqual(meta["n_variants"], 3)

    def test_reference_extracted_from_header(self):
        meta = self._run(self.vcf_path, "chr15_synth")
        self.assertIn("hs37d5", meta["reference"])

    def test_pipeline_and_date_tags_from_filename(self):
        meta = self._run(self.vcf_path, "chr15_synth")
        self.assertEqual(
            meta["pipeline_tag"],
            "phase3_shapeit2_mvncall_integrated_v5b",
        )
        self.assertEqual(meta["date_stamp"], "20130502")

    def test_sample_hash_is_stable_across_runs(self):
        a = self._run(self.vcf_path, "a")
        b = self._run(self.vcf_path, "b")
        self.assertEqual(a["sample_hash"], b["sample_hash"])
        # 32-char md5 hex digest.
        self.assertEqual(len(a["sample_hash"]), 32)

    def test_sample_hash_differs_between_cohorts(self):
        # Build a second VCF whose sample list has one different ID.
        alt_cohort = SyntheticCohort(
            samples=[*default_cohort().samples[:-1], "HG99999"]
        )
        alt_path = standard_filename(self.tmpdir, "22")
        write_vcf(
            alt_path, "22",
            [Variant(pos=16050075, ref="A", alt="G",
                     af_by_pop={"ALL": 0.05})],
            alt_cohort,
        )
        meta_std = self._run(self.vcf_path, "std")
        meta_alt = self._run(alt_path, "alt")
        self.assertNotEqual(meta_std["sample_hash"], meta_alt["sample_hash"])


if __name__ == "__main__":
    unittest.main()
