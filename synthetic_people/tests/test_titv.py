"""Tests for the Ti/Tv calibrator."""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.titv import (  # noqa: E402
    DEFAULT_TARGET_TITV,
    choose_alt,
    is_transition,
    titv_ratio,
)


class TestIsTransition(unittest.TestCase):
    def test_purine_purine(self):
        self.assertTrue(is_transition("A", "G"))
        self.assertTrue(is_transition("G", "A"))

    def test_pyrimidine_pyrimidine(self):
        self.assertTrue(is_transition("C", "T"))
        self.assertTrue(is_transition("T", "C"))

    def test_transversions(self):
        for ref, alt in (("A", "C"), ("A", "T"),
                         ("G", "C"), ("G", "T"),
                         ("C", "A"), ("C", "G"),
                         ("T", "A"), ("T", "G")):
            self.assertFalse(is_transition(ref, alt), (ref, alt))

    def test_case_insensitive(self):
        self.assertTrue(is_transition("a", "g"))
        self.assertTrue(is_transition("c", "T"))


class TestTitvRatio(unittest.TestCase):
    def test_known_counts(self):
        # 4 transitions, 2 transversions → ratio 2.0.
        pairs = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"),
                 ("A", "T"), ("C", "G")]
        self.assertAlmostEqual(titv_ratio(pairs), 2.0)

    def test_skips_indels_and_non_standard(self):
        pairs = [("A", "G"), ("AT", "A"), ("N", "A"),
                 ("C", "T"), ("T", "A")]
        # Only the two real SNVs count: 2 Ti / 1 Tv = 2.0.
        self.assertAlmostEqual(titv_ratio(pairs), 2.0)

    def test_no_transversions(self):
        self.assertEqual(titv_ratio([("A", "G")]), float("inf"))

    def test_empty(self):
        self.assertEqual(titv_ratio([]), 0.0)


class TestChooseAlt(unittest.TestCase):
    def test_never_returns_ref(self):
        rng = random.Random(0)
        for ref in "ACGT":
            for _ in range(200):
                alt = choose_alt(ref, rng)
                self.assertNotEqual(alt, ref)
                self.assertIn(alt, "ACGT")

    def test_rejects_non_standard_ref(self):
        rng = random.Random(0)
        self.assertIsNone(choose_alt("N", rng))
        self.assertIsNone(choose_alt("AT", rng))
        self.assertIsNone(choose_alt("", rng))

    def test_hits_target_ratio_2_1(self):
        """Long-run Ti/Tv should converge on the target."""
        rng = random.Random(42)
        target = 2.1
        pairs = []
        n = 40_000
        for _ in range(n):
            ref = rng.choice("ACGT")
            alt = choose_alt(ref, rng, target=target)
            pairs.append((ref, alt))
        ratio = titv_ratio(pairs)
        # ±5% on the target.
        self.assertAlmostEqual(ratio, target, delta=0.11)

    def test_hits_alternative_targets(self):
        """The calibrator should track whatever target we choose."""
        rng = random.Random(7)
        for target in (0.5, 1.0, 3.0):
            pairs = []
            for _ in range(30_000):
                ref = rng.choice("ACGT")
                pairs.append((ref, choose_alt(ref, rng, target=target)))
            ratio = titv_ratio(pairs)
            self.assertAlmostEqual(ratio, target, delta=max(0.08, 0.05 * target))

    def test_default_target_is_2_1(self):
        self.assertEqual(DEFAULT_TARGET_TITV, 2.1)

    def test_rejects_nonpositive_target(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            choose_alt("A", rng, target=0)
        with self.assertRaises(ValueError):
            choose_alt("A", rng, target=-1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
