"""Distribution checks for syntheticgen.quality.

Runnable either via pytest or bare `python -m unittest`. No numpy
dependency — uses statistics from the stdlib.
"""

from __future__ import annotations

import random
import statistics
import sys
import unittest
from pathlib import Path

# Make the sibling package importable when tests are run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.quality import (  # noqa: E402
    HET_ALT_FRAC,
    ad_from_gt,
    draw_site_quality,
    gq_from_ad,
    poisson,
    sample_lambda,
)


class TestPoisson(unittest.TestCase):
    def test_mean_near_lambda(self):
        rng = random.Random(12345)
        samples = [poisson(30.0, rng) for _ in range(5000)]
        mean = statistics.fmean(samples)
        self.assertAlmostEqual(mean, 30.0, delta=1.0)
        # Variance of Poisson ≈ λ; allow some slack.
        var = statistics.pvariance(samples)
        self.assertAlmostEqual(var, 30.0, delta=3.0)

    def test_floors_at_zero(self):
        rng = random.Random(0)
        samples = [poisson(1.0, rng) for _ in range(2000)]
        self.assertTrue(all(s >= 0 for s in samples))


class TestAD(unittest.TestCase):
    def test_homref_all_ref(self):
        rng = random.Random(1)
        for _ in range(50):
            ref, alt = ad_from_gt("0|0", 30, rng)
            self.assertEqual(ref, 30)
            self.assertEqual(alt, 0)

    def test_homalt_all_alt(self):
        rng = random.Random(2)
        for _ in range(50):
            ref, alt = ad_from_gt("1|1", 40, rng)
            self.assertEqual(ref, 0)
            self.assertEqual(alt, 40)

    def test_het_close_to_half_with_ref_bias(self):
        rng = random.Random(3)
        n = 2000
        alt_fracs = []
        for _ in range(n):
            ref, alt = ad_from_gt("0|1", 30, rng)
            self.assertEqual(ref + alt, 30)
            alt_fracs.append(alt / 30)
        mean_alt_frac = statistics.fmean(alt_fracs)
        # Expect centered around HET_ALT_FRAC with modest tolerance.
        self.assertAlmostEqual(mean_alt_frac, HET_ALT_FRAC, delta=0.03)
        # And it should be visibly below 0.5 (the "perfect" het ratio):
        # reference reads are over-represented on real hets.
        self.assertLess(mean_alt_frac, 0.5)

    def test_ad_sum_equals_dp(self):
        rng = random.Random(4)
        for gt in ("0|0", "0|1", "1|0", "1|1"):
            for dp in (0, 1, 10, 30, 100):
                ref, alt = ad_from_gt(gt, dp, rng)
                self.assertEqual(ref + alt, dp, (gt, dp))


class TestGQ(unittest.TestCase):
    def test_range_clamped(self):
        # High-support cases should go high, contradictions should go low.
        self.assertGreater(gq_from_ad("1|1", 0, 40), 80)
        self.assertGreater(gq_from_ad("0|0", 40, 0), 80)
        # Hom-alt call with zero alt reads → disagreement → low GQ.
        self.assertLess(gq_from_ad("1|1", 40, 0), 10)
        # Depth=0 → GQ=0.
        self.assertEqual(gq_from_ad("0|1", 0, 0), 0)

    def test_het_peaks_near_half(self):
        # GQ for a het call should be highest when AD is roughly 50/50
        # and lower at the tails.
        center = gq_from_ad("0|1", 20, 20)
        skewed = gq_from_ad("0|1", 35, 5)
        self.assertGreater(center, skewed)

    def test_bounds(self):
        for ref in range(0, 30, 5):
            for alt in range(0, 30, 5):
                for gt in ("0|0", "0|1", "1|1"):
                    gq = gq_from_ad(gt, ref, alt)
                    self.assertGreaterEqual(gq, 0)
                    self.assertLessEqual(gq, 99)


class TestDrawSiteQuality(unittest.TestCase):
    def test_tuple_consistency(self):
        rng = random.Random(7)
        for _ in range(200):
            for gt in ("0|0", "0|1", "1|1"):
                dp, ref, alt, gq = draw_site_quality(gt, 30.0, rng)
                self.assertEqual(dp, ref + alt)
                self.assertGreaterEqual(gq, 0)
                self.assertLessEqual(gq, 99)

    def test_cohort_dp_mean(self):
        rng = random.Random(11)
        dps = []
        for _ in range(500):
            lam = sample_lambda(30.0, 3.0, rng)
            for _ in range(10):
                dp, _, _, _ = draw_site_quality("0|1", lam, rng)
                dps.append(dp)
        mean_dp = statistics.fmean(dps)
        self.assertAlmostEqual(mean_dp, 30.0, delta=1.5)

    def test_sample_lambda_clamp(self):
        rng = random.Random(13)
        for _ in range(100):
            lam = sample_lambda(30.0, 3.0, rng)
            self.assertGreaterEqual(lam, 5.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
