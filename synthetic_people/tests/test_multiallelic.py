"""Tests for multi-allelic genotype drawing in background.py."""

from __future__ import annotations

import random
import statistics
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.background import (  # noqa: E402
    alt_dosage,
    alt_dosages,
    phased_gt_from_afs,
)


class TestPhasedGtFromAfs(unittest.TestCase):
    def test_biallelic_reduces_to_bernoulli(self):
        rng = random.Random(42)
        af = 0.3
        counts = {"0|0": 0, "0|1": 0, "1|0": 0, "1|1": 0}
        n = 10000
        for _ in range(n):
            gt = phased_gt_from_afs([af], rng)
            counts[gt] = counts.get(gt, 0) + 1
        # HWE expectations: (1-af)^2 hom-ref, 2*af*(1-af) het, af^2 hom-alt.
        exp_homref = (1 - af) ** 2 * n
        exp_homalt = af ** 2 * n
        homref = counts["0|0"]
        homalt = counts["1|1"]
        self.assertAlmostEqual(homref / n, exp_homref / n, delta=0.02)
        self.assertAlmostEqual(homalt / n, exp_homalt / n, delta=0.02)

    def test_multiallelic_produces_all_alleles(self):
        rng = random.Random(7)
        afs = [0.2, 0.2]  # 40% non-ref total → ~16% het-het, ~4% hom-alt1 etc.
        n = 5000
        # Count occurrences of each allele index on both haplotypes.
        allele_counts = [0, 0, 0]
        for _ in range(n):
            gt = phased_gt_from_afs(afs, rng)
            for tok in gt.split("|"):
                allele_counts[int(tok)] += 1
        total = sum(allele_counts)
        self.assertGreater(allele_counts[1], 0)
        self.assertGreater(allele_counts[2], 0)
        # REF frequency ≈ 0.6.
        self.assertAlmostEqual(allele_counts[0] / total, 0.6, delta=0.02)
        # alt1 + alt2 ≈ 0.4, roughly equally split.
        self.assertAlmostEqual(allele_counts[1] / total, 0.2, delta=0.02)
        self.assertAlmostEqual(allele_counts[2] / total, 0.2, delta=0.02)

    def test_1_2_het_appears_at_realistic_rate(self):
        """With two common alts, 1|2-style hets should show up."""
        rng = random.Random(11)
        afs = [0.3, 0.3]
        het_1_2 = 0
        n = 10000
        for _ in range(n):
            gt = phased_gt_from_afs(afs, rng)
            if gt in ("1|2", "2|1"):
                het_1_2 += 1
        # Expected rate: 2 * af1 * af2 = 0.18.
        self.assertAlmostEqual(het_1_2 / n, 0.18, delta=0.02)


class TestDosages(unittest.TestCase):
    def test_alt_dosage_handles_multi_indices(self):
        self.assertEqual(alt_dosage("0|0"), 0)
        self.assertEqual(alt_dosage("0|1"), 1)
        self.assertEqual(alt_dosage("1|1"), 2)
        self.assertEqual(alt_dosage("0|2"), 1)
        self.assertEqual(alt_dosage("1|2"), 2)
        self.assertEqual(alt_dosage("2|2"), 2)

    def test_per_alt_dosages(self):
        # Biallelic.
        self.assertEqual(alt_dosages("0|1", 1), [1])
        self.assertEqual(alt_dosages("1|1", 1), [2])
        # Multi-allelic.
        self.assertEqual(alt_dosages("0|2", 2), [0, 1])
        self.assertEqual(alt_dosages("1|2", 2), [1, 1])
        self.assertEqual(alt_dosages("2|2", 2), [0, 2])
        # Out-of-range indices are silently clamped (shouldn't happen but
        # the function mustn't crash).
        self.assertEqual(alt_dosages("0|9", 2), [0, 0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
