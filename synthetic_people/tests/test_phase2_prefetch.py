"""Tests for the Phase 2 overlay-loader prefetch.

The CLI submits ClinVar / rsID / COSMIC loaders to a thread pool
*before* the coalescent simulation runs, then resolves the futures
inside the existing overlay block. These tests verify:

* Skip rules — ``--rsid-density 0`` skips the rsID future,
  ``--somatic`` is required for the COSMIC future, default invocation
  schedules ClinVar + rsID and skips COSMIC.
* Identity — the futures resolve to the same data the loader
  functions would have returned if called directly. No serialisation
  / argument-passing surprises.

The real loader functions are bcftools subprocess + I/O against a
gigabyte-scale VCF and would be inappropriate for unit tests, so we
monkey-patch the loader entry points with fakes that record their
arguments and return canned data.
"""

from __future__ import annotations

import sys
import types
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen import cli as cli_module
from syntheticgen.cli import submit_overlays


def _make_args(**overrides) -> types.SimpleNamespace:
    """Build the minimal argparse-style Namespace submit_overlays needs."""
    defaults = dict(
        rsid_density=0.20,
        dbsnp_vcf=None,
        somatic=False,
        cosmic_vcf=None,
    )
    defaults.update(overrides)
    return types.SimpleNamespace(**defaults)


class _FakeLoaders:
    """Captures every call so the tests can assert what was scheduled."""

    def __init__(self):
        self.calls: dict = {
            "clinvar_index": [],
            "rsid_pool": [],
            "cosmic_pool": [],
        }

    def install(self, monkeypatch_target):
        # Replace the loaders the CLI module imported. We patch on the
        # cli module so submit_overlays sees the fakes.
        monkeypatch_target.load_clinvar_index = self._load_clinvar_index
        monkeypatch_target.load_rsid_pool = self._load_rsid_pool
        monkeypatch_target.load_cosmic_records = self._load_cosmic_records

    def _load_clinvar_index(self, vcf, chromosomes, sig_filter,
                            max_per_chrom):
        self.calls["clinvar_index"].append(
            (vcf, tuple(chromosomes), frozenset(sig_filter),
             max_per_chrom))
        return {"CLINVAR_INDEX_FOR": tuple(chromosomes)}

    def _load_rsid_pool(self, vcf, chromosomes, max_per_chrom):
        self.calls["rsid_pool"].append(
            (vcf, tuple(chromosomes), max_per_chrom))
        return {"RSID_POOL_FOR": tuple(chromosomes)}

    def _load_cosmic_records(self, vcf, chromosomes, max_per_chrom):
        self.calls["cosmic_pool"].append(
            (vcf, tuple(chromosomes), max_per_chrom))
        return {"COSMIC_POOL_FOR": tuple(chromosomes)}


class TestSubmitOverlaysSkipRules(unittest.TestCase):
    def setUp(self):
        self.fakes = _FakeLoaders()
        # Stash originals so tearDown restores them — this avoids
        # bleeding fakes into other test modules.
        self._originals = {
            "load_clinvar_index": cli_module.load_clinvar_index,
            "load_rsid_pool": cli_module.load_rsid_pool,
            "load_cosmic_records": cli_module.load_cosmic_records,
        }
        self.fakes.install(cli_module)

    def tearDown(self):
        for name, fn in self._originals.items():
            setattr(cli_module, name, fn)

    def _run(self, args) -> dict:
        with ThreadPoolExecutor(max_workers=3) as ex:
            futures = submit_overlays(
                args, ["21", "22"], Path("/dev/null"),
                {"Pathogenic"}, ex,
            )
        return futures

    def test_default_schedules_clinvar_and_rsid_only(self):
        futures = self._run(_make_args())
        self.assertIsNotNone(futures["clinvar_index"])
        self.assertIsNotNone(futures["rsid_pool"])
        self.assertIsNone(futures["cosmic_pool"])

    def test_rsid_density_zero_skips_rsid(self):
        futures = self._run(_make_args(rsid_density=0))
        self.assertIsNotNone(futures["clinvar_index"])
        self.assertIsNone(futures["rsid_pool"])
        self.assertIsNone(futures["cosmic_pool"])

    def test_somatic_schedules_cosmic(self):
        futures = self._run(_make_args(
            somatic=True, cosmic_vcf=Path("/dev/null")))
        self.assertIsNotNone(futures["cosmic_pool"])

    def test_dbsnp_vcf_used_for_rsid_when_provided(self):
        custom = Path("/tmp/fake-dbsnp.vcf.gz")
        self._run(_make_args(dbsnp_vcf=custom))
        rsid_calls = self.fakes.calls["rsid_pool"]
        self.assertEqual(len(rsid_calls), 1)
        self.assertEqual(rsid_calls[0][0], custom)

    def test_clinvar_used_for_rsid_when_dbsnp_missing(self):
        clinvar_path = Path("/tmp/fake-clinvar.vcf.gz")
        with ThreadPoolExecutor(max_workers=3) as ex:
            submit_overlays(_make_args(), ["22"], clinvar_path,
                            {"Pathogenic"}, ex)
        self.assertEqual(
            self.fakes.calls["rsid_pool"][0][0], clinvar_path)


class TestSubmitOverlaysIdentity(unittest.TestCase):
    """Resolved-future values must match what serial loader calls return."""

    def setUp(self):
        self.fakes = _FakeLoaders()
        self._originals = {
            "load_clinvar_index": cli_module.load_clinvar_index,
            "load_rsid_pool": cli_module.load_rsid_pool,
            "load_cosmic_records": cli_module.load_cosmic_records,
        }
        self.fakes.install(cli_module)

    def tearDown(self):
        for name, fn in self._originals.items():
            setattr(cli_module, name, fn)

    def test_resolved_payload_matches_serial(self):
        args = _make_args(
            somatic=True, cosmic_vcf=Path("/tmp/cosmic.vcf.gz"))
        with ThreadPoolExecutor(max_workers=3) as ex:
            futures = submit_overlays(
                args, ["21", "22"], Path("/tmp/clinvar.vcf.gz"),
                {"Pathogenic"}, ex,
            )
            resolved = {k: (f.result() if f is not None else None)
                        for k, f in futures.items()}

        # The fake loaders return canned dicts keyed on the chromosomes
        # passed in. Verify they match what calling each loader
        # directly would produce.
        self.assertEqual(
            resolved["clinvar_index"],
            cli_module.load_clinvar_index(
                Path("/tmp/clinvar.vcf.gz"), ["21", "22"],
                sig_filter={"Pathogenic"}, max_per_chrom=20_000,
            ),
        )
        self.assertEqual(
            resolved["rsid_pool"],
            cli_module.load_rsid_pool(
                Path("/tmp/clinvar.vcf.gz"), ["21", "22"],
                max_per_chrom=20_000,
            ),
        )
        self.assertEqual(
            resolved["cosmic_pool"],
            cli_module.load_cosmic_records(
                Path("/tmp/cosmic.vcf.gz"), ["21", "22"],
                max_per_chrom=20_000,
            ),
        )

    def test_loader_args_carry_through(self):
        # max_per_chrom is hardcoded inside submit_overlays; sig_filter
        # and chromosomes come from the caller. Pin those.
        with ThreadPoolExecutor(max_workers=3) as ex:
            submit_overlays(
                _make_args(rsid_density=0.5),
                ["19", "20", "21", "22"],
                Path("/tmp/x.vcf.gz"),
                {"Pathogenic", "Likely_pathogenic"},
                ex,
            )

        clinvar_calls = self.fakes.calls["clinvar_index"]
        self.assertEqual(len(clinvar_calls), 1)
        vcf, chroms, sigs, mpc = clinvar_calls[0]
        self.assertEqual(vcf, Path("/tmp/x.vcf.gz"))
        self.assertEqual(chroms, ("19", "20", "21", "22"))
        self.assertEqual(sigs, frozenset(
            {"Pathogenic", "Likely_pathogenic"}))
        self.assertEqual(mpc, 20_000)

        rsid_calls = self.fakes.calls["rsid_pool"]
        self.assertEqual(len(rsid_calls), 1)
        self.assertEqual(rsid_calls[0][1],
                         ("19", "20", "21", "22"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
