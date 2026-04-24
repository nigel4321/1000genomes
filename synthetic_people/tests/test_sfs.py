"""Tests for the SFS sampler."""

from __future__ import annotations

import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.sfs import (  # noqa: E402
    DEFAULT_SFS_ALPHA,
    draw_allele_counts,
    draw_minor_count,
    sfs_histogram,
    singleton_fraction,
    write_sfs_tsv,
)


class TestDrawMinorCount(unittest.TestCase):
    def test_range(self):
        rng = random.Random(0)
        for _ in range(500):
            k = draw_minor_count(40, 2.0, rng)
            self.assertGreaterEqual(k, 1)
            self.assertLessEqual(k, 39)

    def test_rejects_small_n(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            draw_minor_count(1, 2.0, rng)

    def test_rejects_nonpositive_alpha(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            draw_minor_count(10, 0.0, rng)
        with self.assertRaises(ValueError):
            draw_minor_count(10, -1.0, rng)

    def test_alpha_zero_ish_uniformish(self):
        """Very small α should approach uniform over k."""
        rng = random.Random(11)
        counts = {}
        for _ in range(20_000):
            k = draw_minor_count(10, 0.01, rng)
            counts[k] = counts.get(k, 0) + 1
        # All 9 buckets (k=1..9) should see a healthy share.
        for k in range(1, 10):
            self.assertGreater(counts.get(k, 0), 1_000)

    def test_singleton_fraction_default_alpha(self):
        """At default α=2.0 singletons should dominate."""
        rng = random.Random(42)
        n_hap = 40
        draws = [draw_minor_count(n_hap, DEFAULT_SFS_ALPHA, rng)
                 for _ in range(20_000)]
        singletons = sum(1 for k in draws if k == 1)
        self.assertGreater(singletons / len(draws), 0.55)


class TestDrawAlleleCounts(unittest.TestCase):
    def test_biallelic_total_bounded(self):
        rng = random.Random(0)
        for _ in range(200):
            counts = draw_allele_counts(20, 1, 2.0, rng)
            self.assertEqual(len(counts), 1)
            self.assertGreaterEqual(counts[0], 1)
            self.assertLessEqual(counts[0], 19)

    def test_multiallelic_total_bounded(self):
        rng = random.Random(1)
        for _ in range(200):
            counts = draw_allele_counts(20, 3, 2.0, rng)
            self.assertEqual(len(counts), 3)
            self.assertLessEqual(sum(counts), 19)
            for c in counts:
                self.assertGreaterEqual(c, 1)

    def test_rejects_invalid_n_alts(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            draw_allele_counts(10, 0, 2.0, rng)


class TestSfsHistogram(unittest.TestCase):
    def test_collects_per_alt(self):
        sites = [
            {"acs": [3, 1]},
            {"acs": [1]},
            {"acs": [5]},
            {"acs": [1, 1]},
        ]
        hist = sfs_histogram(sites)
        self.assertEqual(hist, {1: 4, 3: 1, 5: 1})

    def test_singleton_fraction(self):
        hist = {1: 60, 2: 30, 3: 10}
        self.assertAlmostEqual(singleton_fraction(hist), 0.6)

    def test_singleton_fraction_empty(self):
        self.assertEqual(singleton_fraction({}), 0.0)

    def test_write_tsv_round_trip(self):
        hist = {1: 100, 2: 40, 5: 8}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "summary" / "sfs.tsv"
            write_sfs_tsv(path, hist)
            content = path.read_text().strip().splitlines()
        self.assertEqual(content[0], "ac\tn_sites")
        self.assertIn("1\t100", content)
        self.assertIn("2\t40", content)
        self.assertIn("5\t8", content)


if __name__ == "__main__":
    unittest.main(verbosity=2)
