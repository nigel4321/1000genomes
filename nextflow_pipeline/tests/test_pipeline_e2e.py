"""End-to-end: invoke `nextflow run main.nf` against synthetic VCFs.

Generates a small chr15 + chr22 cohort, runs the pipeline, and verifies the
three published reports + carriers.tsv look right. Skipped automatically if
nextflow is not on PATH.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    PIPELINE_DIR,
    default_cohort,
    in_range_variants,
    require_tools,
    standard_filename,
)
from synthetic_vcf import Variant, write_vcf


class PipelineE2ETest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        require_tools("nextflow", "bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="pipeline_e2e_")
        cls.input_dir = os.path.join(cls.tmpdir, "vcfs")
        os.makedirs(cls.input_dir)

        # chr15: target variant in range.
        write_vcf(
            standard_filename(cls.input_dir, "15"),
            "15", in_range_variants(), default_cohort(),
        )
        # chr22: does not contain chr15 → should classify as not_applicable.
        write_vcf(
            standard_filename(cls.input_dir, "22"),
            "22",
            [Variant(pos=16050075, ref="A", alt="G",
                     af_by_pop={"ALL": 0.05})],
            default_cohort(),
        )

        cls.outdir = os.path.join(cls.tmpdir, "results")
        cls.workdir = os.path.join(cls.tmpdir, "work")

    @classmethod
    def tearDownClass(cls):
        # Leave the work dir behind if NF_KEEP is set, so a failing run can
        # be inspected. Otherwise clean up — it is tens of MB per run.
        if not os.environ.get("NF_KEEP"):
            shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_pipeline_runs_and_produces_reports(self):
        glob = os.path.join(self.input_dir, "*.vcf.gz")
        cmd = [
            "nextflow", "run", "main.nf",
            "--input", glob,
            "--outdir", self.outdir,
            "-work-dir", self.workdir,
            "--variant_name", "rs12913832",
            "--variant_chrom", "15",
            "--variant_pos", "28365618",
            "--variant_ref", "A",
            "--variant_alt", "G",
            "--variant_min_af", "0.05",
            "--variant_max_af", "1.0",
            "-ansi-log", "false",
        ]
        proc = subprocess.run(cmd, cwd=PIPELINE_DIR,
                              capture_output=True, text=True)
        self.assertEqual(
            proc.returncode, 0,
            f"nextflow failed\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}",
        )

        # All three reports + carriers.tsv should be published.
        for name in ("metadata_report.md", "variant_report.md",
                     "carriers_report.md", "carriers.tsv"):
            self.assertTrue(
                os.path.isfile(os.path.join(self.outdir, name)),
                f"missing output: {name}",
            )

        with open(os.path.join(self.outdir, "variant_report.md")) as fh:
            variant_md = fh.read()
        # chr15 file should be in-range, chr22 file should be not_applicable.
        self.assertIn("present_in_range", variant_md)
        self.assertIn("not_applicable", variant_md)
        self.assertIn("rs12913832", variant_md)

        with open(os.path.join(self.outdir, "carriers_report.md")) as fh:
            carrier_md = fh.read()
        self.assertIn("Allele-count integrity check", carrier_md)
        self.assertIn("Heterozygotes", carrier_md)

        with open(os.path.join(self.outdir, "carriers.tsv")) as fh:
            carriers_tsv = fh.read().splitlines()
        header = carriers_tsv[0].split("\t")
        self.assertEqual(header[:3], ["file", "sample", "variant_id"])
        data_rows = [l for l in carriers_tsv[1:] if l.strip()]
        # All 4 EUR samples should be homozygous alt → at minimum 4 carriers.
        self.assertGreaterEqual(len(data_rows), 4)
        samples_in_tsv = {l.split("\t")[1] for l in data_rows}
        self.assertTrue(
            {"HG00096", "HG00097", "HG00099", "HG00100"}.issubset(samples_in_tsv),
        )

        with open(os.path.join(self.outdir, "metadata_report.md")) as fh:
            meta_md = fh.read()
        self.assertIn("**Files scanned:** 2", meta_md)
        self.assertIn("phase3_shapeit2_mvncall_integrated_v5b", meta_md)


if __name__ == "__main__":
    unittest.main()
