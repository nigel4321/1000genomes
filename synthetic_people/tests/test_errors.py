"""Tests for the lightweight sequencing-error model (M9).

Pure-Python — no msprime / bcftools dependency, so this file runs in
any environment alongside the M7/M8 tests.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.errors import (
    DEFAULT_DROPOUT_RATE,
    DEFAULT_GT_ERROR_RATE,
    maybe_dropout,
    maybe_flip_gt,
    merge_stats,
    new_error_stats,
)


class TestMaybeFlipGt(unittest.TestCase):
    def test_zero_rate_is_noop(self):
        rng = random.Random(0)
        for gt in ("0|0", "0|1", "1|0", "1|1", "1|2"):
            new_gt, flipped = maybe_flip_gt(gt, rng, 0.0)
            self.assertEqual(new_gt, gt)
            self.assertFalse(flipped)

    def test_negative_rate_is_noop(self):
        rng = random.Random(0)
        new_gt, flipped = maybe_flip_gt("0|1", rng, -0.5)
        self.assertEqual(new_gt, "0|1")
        self.assertFalse(flipped)

    def test_full_rate_always_flips(self):
        """At rate=1.0 every call should flip to a different GT."""
        rng = random.Random(1)
        for gt in ("0|0", "0|1", "1|0", "1|1"):
            for _ in range(20):
                new_gt, flipped = maybe_flip_gt(gt, rng, 1.0)
                self.assertTrue(flipped)
                self.assertNotEqual(new_gt, gt)

    def test_realised_rate_close_to_target(self):
        """A 1% rate over 10k draws should land in [~50, ~150] flips
        (3-sigma binomial band around 100)."""
        rng = random.Random(42)
        n = 10_000
        flipped_count = sum(
            1 for _ in range(n)
            if maybe_flip_gt("0|0", rng, 0.01)[1]
        )
        self.assertGreater(flipped_count, 50)
        self.assertLess(flipped_count, 150)

    def test_biallelic_flip_targets_are_valid(self):
        """Every flip from a biallelic GT must land on another
        biallelic GT in the canonical set."""
        rng = random.Random(7)
        valid = {"0|0", "0|1", "1|0", "1|1"}
        for gt in ("0|0", "0|1", "1|0", "1|1"):
            for _ in range(200):
                new_gt, flipped = maybe_flip_gt(gt, rng, 1.0)
                self.assertIn(new_gt, valid)
                self.assertTrue(flipped)
                self.assertNotEqual(new_gt, gt)

    def test_hom_to_het_dominates(self):
        """Hom-ref → het should be the most common flip (0.7 weight)."""
        rng = random.Random(2)
        n = 2000
        targets = {"0|1": 0, "1|1": 0}
        for _ in range(n):
            new_gt, _ = maybe_flip_gt("0|0", rng, 1.0)
            if new_gt in targets:
                targets[new_gt] += 1
        # Wide tolerance: 0.7±0.05 → 1300–1500
        self.assertGreater(targets["0|1"], 1200)
        self.assertLess(targets["0|1"], 1600)
        self.assertLess(targets["1|1"], targets["0|1"])

    def test_het_splits_evenly(self):
        """0|1 → 0|0 vs 1|1 should be ~50/50."""
        rng = random.Random(3)
        n = 2000
        a = b = 0
        for _ in range(n):
            new_gt, _ = maybe_flip_gt("0|1", rng, 1.0)
            if new_gt == "0|0":
                a += 1
            elif new_gt == "1|1":
                b += 1
        # 50/50 ± 5%
        self.assertGreater(a, 800)
        self.assertGreater(b, 800)

    def test_multiallelic_collapses_to_ref(self):
        """`1|2`-style hets should always carry a 0 on one side after
        flipping (caller drops one allele)."""
        rng = random.Random(4)
        for _ in range(100):
            new_gt, flipped = maybe_flip_gt("1|2", rng, 1.0)
            self.assertTrue(flipped)
            self.assertIn("0", new_gt.split("|"))

    def test_unparseable_gt_passes_through(self):
        rng = random.Random(0)
        # Even at rate=1.0 a malformed token shouldn't blow up.
        new_gt, flipped = maybe_flip_gt("./.", rng, 1.0)
        self.assertEqual(new_gt, "./.")
        self.assertFalse(flipped)
        new_gt, flipped = maybe_flip_gt("0", rng, 1.0)
        self.assertEqual(new_gt, "0")
        self.assertFalse(flipped)

    def test_seed_reproducibility(self):
        rng_a = random.Random(99)
        rng_b = random.Random(99)
        gts = ("0|0", "0|1", "1|0", "1|1")
        for gt in gts * 50:
            self.assertEqual(
                maybe_flip_gt(gt, rng_a, 0.05),
                maybe_flip_gt(gt, rng_b, 0.05),
            )


class TestMaybeDropout(unittest.TestCase):
    def test_zero_rate_never_drops(self):
        rng = random.Random(0)
        for _ in range(1000):
            self.assertFalse(maybe_dropout(rng, 0.0))

    def test_full_rate_always_drops(self):
        rng = random.Random(0)
        for _ in range(100):
            self.assertTrue(maybe_dropout(rng, 1.0))

    def test_realised_rate_close_to_target(self):
        rng = random.Random(123)
        n = 10_000
        drops = sum(1 for _ in range(n) if maybe_dropout(rng, 0.01))
        self.assertGreater(drops, 50)
        self.assertLess(drops, 150)

    def test_seed_reproducibility(self):
        a = random.Random(5)
        b = random.Random(5)
        for _ in range(500):
            self.assertEqual(maybe_dropout(a, 0.1),
                             maybe_dropout(b, 0.1))


class TestStats(unittest.TestCase):
    def test_new_error_stats_has_expected_keys(self):
        s = new_error_stats()
        self.assertEqual(s, {"flipped": 0, "dropped": 0,
                             "total_calls": 0})

    def test_merge_stats_adds_in_place(self):
        a = new_error_stats()
        a["flipped"] = 3
        a["total_calls"] = 100
        b = {"flipped": 2, "dropped": 1, "total_calls": 50}
        merge_stats(a, b)
        self.assertEqual(a, {"flipped": 5, "dropped": 1,
                             "total_calls": 150})

    def test_merge_stats_seeds_missing_keys(self):
        a = {}
        merge_stats(a, {"x": 7})
        self.assertEqual(a, {"x": 7})


class TestDefaultConstants(unittest.TestCase):
    def test_defaults_lock_in(self):
        # Lock the documented defaults so unintentional changes surface
        # as test failures.
        self.assertEqual(DEFAULT_GT_ERROR_RATE, 0.001)
        self.assertEqual(DEFAULT_DROPOUT_RATE, 0.0005)


if __name__ == "__main__":
    unittest.main(verbosity=2)
