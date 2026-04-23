"""Tests for the three markdown-report builders under bin/.

These tests bypass Nextflow: they generate inputs (meta.json /
variant.json / carriers.tsv) into a temp dir, then invoke each
build_*_report.py script and inspect the emitted markdown / TSV.
"""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import bin_script, run_script


class MetadataReportTest(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="metadata_report_test_")
        self.input_dir = os.path.join(self.tmpdir, "_metas")
        os.makedirs(self.input_dir)

    def _write_meta(self, name: str, **overrides):
        meta = {
            "name": name,
            "file": f"{name}.vcf.gz",
            "contigs": ["15"],
            "n_variants": 1000,
            "n_samples": 20,
            "sample_hash": "abc123",
            "reference": "hs37d5",
            "size_bytes": 1024 * 1024,
            "pipeline_tag": "phase3_shapeit2_mvncall_integrated_v5b",
            "date_stamp": "20130502",
        }
        meta.update(overrides)
        path = os.path.join(self.input_dir, f"{name}.meta.json")
        with open(path, "w") as fh:
            json.dump(meta, fh)
        return path

    def _build(self) -> str:
        out = os.path.join(self.tmpdir, "metadata_report.md")
        run_script(
            bin_script("build_metadata_report.py"),
            "--input-dir", self.input_dir,
            "--output", out,
        )
        with open(out) as fh:
            return fh.read()

    def test_single_cohort_report(self):
        self._write_meta("chr15")
        self._write_meta("chr22", contigs=["22"])
        md = self._build()
        self.assertIn("VCF Cohort Metadata Report", md)
        self.assertIn("**Files scanned:** 2", md)
        self.assertIn("**Total variants across files:** 2,000", md)
        self.assertIn("All files share the same 20-sample cohort", md)
        self.assertIn("chr15.vcf.gz", md)
        self.assertIn("chr22.vcf.gz", md)

    def test_cohort_mismatch_flagged(self):
        self._write_meta("chr15", sample_hash="aaa")
        self._write_meta("chr22", sample_hash="bbb", contigs=["22"])
        md = self._build()
        self.assertIn("Sample sets differ", md)

    def test_partial_variant_counts_noted(self):
        # Older tabix indices may lack count metadata → n_variants=None.
        self._write_meta("chr15")
        self._write_meta("chr22", n_variants=None, contigs=["22"])
        md = self._build()
        self.assertIn("partial", md.lower())


class VariantReportTest(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="variant_report_test_")
        self.input_dir = os.path.join(self.tmpdir, "_results")
        os.makedirs(self.input_dir)

    def _write_result(self, name: str, **overrides):
        default = {
            "name": name, "file": f"{name}.vcf.gz",
            "variant_chrom": "15", "variant_pos": 28365618,
            "variant_ref": "A", "variant_alt": "G",
            "min_af": 0.05, "max_af": 1.0,
            "status": "present_in_range",
            "variant_id": "rs12913832",
            "chrom": "15", "pos": 28365618, "ref": "A", "alt": "G",
            "ac": 10, "an": 40, "af": 0.25,
            "eur_af": 1.0, "afr_af": 0.0, "eas_af": 0.125,
            "sas_af": 0.25, "amr_af": 0.0,
            "n_samples_scanned": 20, "n_carriers": 8,
            "n_heterozygotes": 2, "n_homozygous_alt": 4,
        }
        default.update(overrides)
        path = os.path.join(self.input_dir, f"{name}.variant.json")
        with open(path, "w") as fh:
            json.dump(default, fh)
        return path

    def _build(self, min_af=0.05, max_af=1.0) -> str:
        out = os.path.join(self.tmpdir, "variant_report.md")
        run_script(
            bin_script("build_variant_report.py"),
            "--input-dir", self.input_dir,
            "--output", out,
            "--variant-name", "rs12913832",
            "--variant-chrom", "15",
            "--variant-pos", "28365618",
            "--variant-ref", "A",
            "--variant-alt", "G",
            "--min-af", str(min_af),
            "--max-af", str(max_af),
        )
        with open(out) as fh:
            return fh.read()

    def test_lists_in_range_files(self):
        self._write_result("chr15")
        md = self._build()
        self.assertIn("rs12913832", md)
        self.assertIn("15:28365618", md)
        self.assertIn("chr15.vcf.gz", md)
        self.assertIn("present_in_range", md)

    def test_mixed_statuses_summarised(self):
        self._write_result("hit", status="present_in_range")
        self._write_result("miss", status="absent",
                           af=None, ac=None, an=None,
                           note="no A>G at position")
        self._write_result("na", status="not_applicable",
                           af=None, ac=None, an=None,
                           note="chromosome '15' not present in VCF")
        md = self._build()
        self.assertIn("**present_in_range:** 1", md)
        self.assertIn("**absent:** 1", md)
        self.assertIn("**not_applicable:** 1", md)

    def test_empty_in_range_section_when_no_hits(self):
        self._write_result("miss", status="absent",
                           af=None, ac=None, an=None,
                           note="no A>G at position")
        md = self._build()
        self.assertIn("No files contain the variant", md)


class CarrierReportTest(unittest.TestCase):

    HEADER = (
        "file\tsample\tvariant_id\tchrom\tpos\tref\talt\tgenotype\talt_dosage"
    )

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="carrier_report_test_")
        self.input_dir = os.path.join(self.tmpdir, "_carriers")
        os.makedirs(self.input_dir)

    def _write_carrier_file(self, name: str, rows: list[tuple]):
        path = os.path.join(self.input_dir, f"{name}.carriers.tsv")
        with open(path, "w") as fh:
            fh.write(self.HEADER + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        return path

    def _build(self) -> tuple[str, str]:
        tsv = os.path.join(self.tmpdir, "carriers.tsv")
        md = os.path.join(self.tmpdir, "carriers_report.md")
        run_script(
            bin_script("build_carrier_report.py"),
            "--input-dir", self.input_dir,
            "--output-tsv", tsv,
            "--output-md", md,
            "--variant-name", "rs12913832",
        )
        with open(tsv) as fh:
            tsv_text = fh.read()
        with open(md) as fh:
            md_text = fh.read()
        return tsv_text, md_text

    def _row(self, file_name: str, sample: str, gt: str, dosage: int):
        return (file_name, sample, "rs12913832", "15", "28365618",
                "A", "G", gt, dosage)

    def test_aggregates_carriers_and_integrity_check(self):
        # 2 heterozygotes + 3 homozygotes → AC = 2 + 2*3 = 8
        rows = [
            self._row("chr15.vcf.gz", "HG00096", "0|1", 1),
            self._row("chr15.vcf.gz", "HG00097", "1|0", 1),
            self._row("chr15.vcf.gz", "HG00099", "1|1", 2),
            self._row("chr15.vcf.gz", "HG00100", "1|1", 2),
            self._row("chr15.vcf.gz", "HG00403", "1|1", 2),
        ]
        self._write_carrier_file("chr15", rows)
        tsv, md = self._build()
        # TSV must contain exactly the 5 carrier rows + header.
        self.assertEqual(len(tsv.strip().splitlines()), 6)
        self.assertIn("**Heterozygotes (dosage = 1):** 2", md)
        self.assertIn("**Homozygous alt (dosage = 2):** 3", md)
        self.assertIn("het + 2·hom = 2 + 2·3 = **8**", md)
        self.assertIn("**Unique individuals carrying the alt allele:** 5", md)

    def test_handles_header_only_files_silently(self):
        # This is the shape emitted when a file is not_applicable/absent.
        self._write_carrier_file("chr22", [])
        self._write_carrier_file("chr15", [
            self._row("chr15.vcf.gz", "HG00096", "0|1", 1),
        ])
        tsv, md = self._build()
        self.assertEqual(len(tsv.strip().splitlines()), 2)  # header + 1 row
        self.assertIn("**Heterozygotes (dosage = 1):** 1", md)

    def test_per_file_breakdown(self):
        self._write_carrier_file("chr15", [
            self._row("chr15.vcf.gz", "HG00096", "0|1", 1),
            self._row("chr15.vcf.gz", "HG00099", "1|1", 2),
        ])
        self._write_carrier_file("chr22", [
            self._row("chr22.vcf.gz", "HG00403", "0|1", 1),
        ])
        _, md = self._build()
        self.assertIn("chr15.vcf.gz", md)
        self.assertIn("chr22.vcf.gz", md)
        # Per-file row for chr15: 1 het + 1 hom = 2
        self.assertIn("| `chr15.vcf.gz` | 1 | 1 | 2 |", md)
        self.assertIn("| `chr22.vcf.gz` | 1 | 0 | 1 |", md)


class QcReportTest(unittest.TestCase):
    """build_qc_report.py — aggregate qc.json files into qc_report.md."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="qc_report_test_")
        self.input_dir = os.path.join(self.tmpdir, "_qc")
        os.makedirs(self.input_dir)

    def _write(self, name: str, *, pass_: bool = True,
               errors=None, warnings=None, **check_overrides):
        default_checks = {
            "size_bytes": 2048,
            "tabix_index_present": True,
            "header_parseable": True,
            "reference": "hs37d5",
            "reference_build": "GRCh37",
            "sample_count": 20,
            "contigs": ["15"],
            "human_contigs": ["15"],
            "non_human_contigs": [],
            "n_variants": 100,
            "has_format_gt": True,
            "info_fields": ["AC", "AF", "AN"],
            "info_has_af": True,
            "info_has_ac_an": True,
        }
        default_checks.update(check_overrides)
        payload = {
            "name": name,
            "file": f"{name}.vcf.gz",
            "pass": pass_,
            "errors": errors or [],
            "warnings": warnings or [],
            "checks": default_checks,
        }
        path = os.path.join(self.input_dir, f"{name}.qc.json")
        with open(path, "w") as fh:
            json.dump(payload, fh)
        return path

    def _build(self) -> str:
        out = os.path.join(self.tmpdir, "qc_report.md")
        run_script(
            bin_script("build_qc_report.py"),
            "--input-dir", self.input_dir,
            "--output", out,
        )
        with open(out) as fh:
            return fh.read()

    def test_all_pass_summary(self):
        self._write("chr15")
        self._write("chr22", contigs=["22"], human_contigs=["22"])
        md = self._build()
        self.assertIn("**Files scanned:** 2", md)
        self.assertIn("**Passed:** 2", md)
        self.assertIn("**Failed:** 0", md)
        self.assertNotIn("## Hard failures", md)
        self.assertNotIn("## Warnings", md)

    def test_hard_failure_rendered(self):
        self._write("chr15")
        self._write(
            "broken", pass_=False,
            errors=["missing tabix index: expected broken.vcf.gz.tbi"],
            tabix_index_present=False,
        )
        md = self._build()
        self.assertIn("**Failed:** 1", md)
        self.assertIn("## Hard failures", md)
        self.assertIn("`broken.vcf.gz`", md)
        self.assertIn("missing tabix index", md)

    def test_warnings_rendered(self):
        self._write(
            "weird", pass_=True,
            warnings=[
                "reference 'mm10' does not match any recognised human build",
                "no indexed contigs match human chromosome names",
            ],
            reference="mm10", reference_build=None,
            contigs=["chr_mouse"], human_contigs=[],
            non_human_contigs=["chr_mouse"],
        )
        md = self._build()
        self.assertIn("## Warnings", md)
        self.assertIn("`weird.vcf.gz`", md)
        self.assertIn("mm10", md)

    def test_per_file_detail_table(self):
        self._write("chr15")
        md = self._build()
        self.assertIn("| File | Status | Samples | Variants | Reference", md)
        self.assertIn("| `chr15.vcf.gz` | PASS | 20 | 100 ", md)


if __name__ == "__main__":
    unittest.main()
