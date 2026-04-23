"""Tests for bin/qc_validate.py — pass, warn, and hard-fail paths."""

import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    bin_script,
    default_cohort,
    in_range_variants,
    require_tools,
    standard_filename,
)
from synthetic_vcf import Variant, write_vcf


def _run_qc(vcf_path: str, name: str, out_dir: str, strict: bool = False):
    out_json = os.path.join(out_dir, f"{name}.qc.json")
    cmd = [
        sys.executable, bin_script("qc_validate.py"),
        "--vcf", vcf_path, "--name", name, "--out", out_json,
    ]
    if strict:
        cmd.append("--strict")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    result = None
    if os.path.isfile(out_json):
        with open(out_json) as fh:
            result = json.load(fh)
    return proc, result


class QcValidatePassTest(unittest.TestCase):
    """A valid human-VCF should pass with no errors or warnings."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="qc_pass_")
        cls.vcf = write_vcf(
            standard_filename(cls.tmpdir, "15"), "15",
            in_range_variants(), default_cohort(),
        )

    def test_clean_file_passes(self):
        proc, result = _run_qc(self.vcf, "chr15", self.tmpdir, strict=True)
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        self.assertTrue(result["pass"])
        self.assertEqual(result["errors"], [])
        self.assertEqual(result["warnings"], [])

    def test_checks_capture_human_attributes(self):
        _, result = _run_qc(self.vcf, "chr15", self.tmpdir)
        c = result["checks"]
        self.assertEqual(c["reference_build"], "GRCh37")
        self.assertEqual(c["human_contigs"], ["15"])
        self.assertEqual(c["non_human_contigs"], [])
        self.assertTrue(c["has_format_gt"])
        self.assertTrue(c["info_has_af"])
        self.assertEqual(c["sample_count"], 20)


class QcValidateWarningsTest(unittest.TestCase):
    """Soft checks should raise warnings without flipping pass=false."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="qc_warn_")

    def _make_vcf(self, chrom="15", **kwargs):
        path = standard_filename(self.tmpdir, chrom)
        write_vcf(
            path, chrom,
            [Variant(pos=28365618, ref="A", alt="G",
                     af_by_pop={"ALL": 0.1})],
            default_cohort(),
            **kwargs,
        )
        return path

    def test_non_human_reference_warns(self):
        vcf = self._make_vcf(reference="mm10.fa")  # mouse reference
        _, result = _run_qc(vcf, "mouse", self.tmpdir)
        self.assertTrue(result["pass"])
        self.assertIsNone(result["checks"]["reference_build"])
        self.assertTrue(
            any("does not match" in w for w in result["warnings"]),
            msg=result["warnings"],
        )

    def test_non_human_contig_warns(self):
        # Emit data on a non-standard contig and override the header to
        # declare only that contig (otherwise bcftools complains about
        # records on an undeclared contig).
        vcf = self._make_vcf(
            chrom="synthetic_scaffold_1",
            contigs_override={"synthetic_scaffold_1": 1_000_000},
        )
        _, result = _run_qc(vcf, "scaffold", self.tmpdir)
        self.assertTrue(result["pass"])
        self.assertEqual(result["checks"]["human_contigs"], [])
        self.assertIn("synthetic_scaffold_1",
                      result["checks"]["non_human_contigs"])
        self.assertTrue(
            any("human chromosome" in w for w in result["warnings"]),
            msg=result["warnings"],
        )

    def test_missing_af_and_acan_warns(self):
        # Drop all three so the pipeline's fill-tags fallback would be the
        # only thing left — QC should flag this up-front.
        vcf = self._make_vcf(info_declarations=[
            '##INFO=<ID=NS,Number=1,Type=Integer,Description="ns">',
            '##INFO=<ID=VT,Number=.,Type=String,Description="vt">',
        ])
        _, result = _run_qc(vcf, "no_af", self.tmpdir)
        self.assertTrue(result["pass"])
        self.assertFalse(result["checks"]["info_has_af"])
        self.assertFalse(result["checks"]["info_has_ac_an"])
        self.assertTrue(
            any("AF" in w and "AC" in w for w in result["warnings"]),
            msg=result["warnings"],
        )

    def test_missing_gt_format_warns(self):
        vcf = self._make_vcf(declare_format_gt=False)
        _, result = _run_qc(vcf, "no_gt", self.tmpdir)
        self.assertTrue(result["pass"])
        self.assertFalse(result["checks"]["has_format_gt"])
        self.assertTrue(
            any("GT" in w for w in result["warnings"]),
            msg=result["warnings"],
        )


class QcValidateLegacyIndexTest(unittest.TestCase):
    """Older tabix indices (e.g. 1000G Phase 3) lack count metadata so
    `bcftools index -s` exits 1 with empty stdout. QC must fall back to
    `tabix -l` — before the fix, real chr20 was flagged as "no contigs"."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="qc_legacy_idx_")
        cls.vcf = write_vcf(
            standard_filename(cls.tmpdir, "15"), "15",
            in_range_variants(), default_cohort(),
        )

    def _shadowed_bcftools(self) -> str:
        """Create a bcftools wrapper that simulates a legacy tabix index.

        For `bcftools index -s`: exit 1, empty stdout (matches what htslib
        emits for pre-counts indices). Everything else delegates to the real
        binary.
        """
        real = os.environ.get("REAL_BCFTOOLS") or \
            subprocess.run(["which", "bcftools"], capture_output=True,
                           text=True).stdout.strip()
        shadow_dir = os.path.join(self.tmpdir, "shadow_bin")
        os.makedirs(shadow_dir, exist_ok=True)
        wrapper = os.path.join(shadow_dir, "bcftools")
        with open(wrapper, "w") as fh:
            fh.write(
                "#!/bin/bash\n"
                f'real={real}\n'
                'if [[ "$1" == "index" && "$2" == "-s" ]]; then\n'
                '  exit 1\n'
                'fi\n'
                'exec "$real" "$@"\n'
            )
        os.chmod(wrapper, 0o755)
        return shadow_dir

    def test_falls_back_to_tabix_l(self):
        shadow = self._shadowed_bcftools()
        env = os.environ.copy()
        env["PATH"] = shadow + os.pathsep + env["PATH"]
        out = os.path.join(self.tmpdir, "legacy.qc.json")
        proc = subprocess.run(
            [sys.executable, bin_script("qc_validate.py"),
             "--vcf", self.vcf, "--name", "legacy",
             "--out", out, "--strict"],
            env=env, capture_output=True, text=True,
        )
        self.assertEqual(
            proc.returncode, 0,
            msg=f"expected pass via tabix -l fallback\n"
                f"stderr:\n{proc.stderr}",
        )
        with open(out) as fh:
            result = json.load(fh)
        self.assertTrue(result["pass"])
        self.assertEqual(result["checks"]["contigs"], ["15"])
        # Counts unavailable from the shadow → n_variants should be null.
        self.assertIsNone(result["checks"]["n_variants"])


class QcValidateHardFailuresTest(unittest.TestCase):
    """Hard checks abort in strict mode and are reported otherwise."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="qc_fail_")

    def test_missing_file_fails(self):
        missing = os.path.join(self.tmpdir, "does_not_exist.vcf.gz")
        proc, result = _run_qc(missing, "ghost", self.tmpdir, strict=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(result["pass"])
        self.assertTrue(
            any("does not exist" in e for e in result["errors"]),
            msg=result["errors"],
        )

    def test_missing_tabix_index_fails(self):
        vcf = write_vcf(
            standard_filename(self.tmpdir, "15"), "15",
            in_range_variants(), default_cohort(),
        )
        # Remove the index the generator created.
        os.remove(vcf + ".tbi")
        proc, result = _run_qc(vcf, "no_idx", self.tmpdir, strict=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(result["pass"])
        self.assertTrue(
            any("tabix index" in e for e in result["errors"]),
            msg=result["errors"],
        )

    def test_empty_file_fails(self):
        empty = os.path.join(self.tmpdir, "empty.vcf.gz")
        open(empty, "w").close()
        # Also write an empty .tbi so the missing-index check doesn't fire
        # first — we want this test to isolate the empty-file path.
        open(empty + ".tbi", "w").close()
        proc, result = _run_qc(empty, "empty", self.tmpdir, strict=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(result["pass"])
        self.assertTrue(
            any("empty" in e for e in result["errors"]),
            msg=result["errors"],
        )

    def test_non_strict_does_not_exit_nonzero_even_on_hard_error(self):
        missing = os.path.join(self.tmpdir, "does_not_exist.vcf.gz")
        proc, result = _run_qc(missing, "ghost", self.tmpdir, strict=False)
        self.assertEqual(proc.returncode, 0)
        self.assertFalse(result["pass"])
        self.assertTrue(result["errors"])


if __name__ == "__main__":
    unittest.main()
