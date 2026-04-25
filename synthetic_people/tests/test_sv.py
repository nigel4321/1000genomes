"""Tests for the structural-variant generator (M8).

Pure-Python — no msprime / bcftools dependency, so this file runs in
any environment alongside the M7 overlay tests.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.sv import (
    DEFAULT_SV_LENGTH_MAX_BP,
    DEFAULT_SV_LENGTH_MIN_BP,
    SV_TYPES,
    _build_sv_record,
    _draw_length,
    generate_person_svs,
)


class TestDrawLength(unittest.TestCase):
    def test_within_bounds(self):
        rng = random.Random(0)
        for _ in range(500):
            v = _draw_length(rng, 50, 10_000)
            self.assertGreaterEqual(v, 50)
            self.assertLessEqual(v, 10_000)

    def test_log_uniform_skews_toward_short(self):
        """A log-uniform draw on [50, 10_000] should put more mass on
        short SVs than long. Median should land well below the
        arithmetic midpoint of 5_025 bp."""
        rng = random.Random(1)
        draws = [_draw_length(rng, 50, 10_000) for _ in range(2_000)]
        draws.sort()
        median = draws[len(draws) // 2]
        # log-midpoint is sqrt(50 * 10_000) ≈ 707; allow a wide window.
        self.assertLess(median, 1_500)

    def test_collapsed_range_returns_endpoint(self):
        rng = random.Random(0)
        self.assertEqual(_draw_length(rng, 100, 100), 100)

    def test_rejects_invalid_range(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            _draw_length(rng, 0, 100)
        with self.assertRaises(ValueError):
            _draw_length(rng, 100, 50)


class TestBuildSvRecord(unittest.TestCase):
    def test_del_has_negative_svlen(self):
        rng = random.Random(0)
        rec = _build_sv_record("22", 1000, 500, "DEL", "0|1", rng)
        self.assertEqual(rec["svtype"], "DEL")
        self.assertEqual(rec["svlen"], -500)
        self.assertEqual(rec["end"], 1500)
        self.assertEqual(rec["alts"], ["<DEL>"])

    def test_dup_has_positive_svlen(self):
        rng = random.Random(0)
        rec = _build_sv_record("22", 1000, 500, "DUP", "1|1", rng)
        self.assertEqual(rec["svtype"], "DUP")
        self.assertEqual(rec["svlen"], 500)
        self.assertEqual(rec["alts"], ["<DUP>"])
        self.assertEqual(rec["gt"], "1|1")

    def test_inv_has_positive_svlen(self):
        rng = random.Random(0)
        rec = _build_sv_record("22", 1000, 750, "INV", "0|1", rng)
        self.assertEqual(rec["svtype"], "INV")
        self.assertEqual(rec["svlen"], 750)
        self.assertEqual(rec["alts"], ["<INV>"])

    def test_anchor_base_is_standard(self):
        rng = random.Random(0)
        rec = _build_sv_record("22", 1000, 100, "DEL", "0|1", rng)
        self.assertIn(rec["ref"], {"A", "C", "G", "T"})
        self.assertEqual(len(rec["ref"]), 1)

    def test_cipos_default_imprecise(self):
        rng = random.Random(0)
        rec = _build_sv_record("22", 1000, 500, "DEL", "0|1", rng)
        self.assertEqual(rec["cipos"], (-50, 50))

    def test_unknown_svtype_raises(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            _build_sv_record("22", 1000, 500, "BND", "0|1", rng)


class TestGeneratePersonSvs(unittest.TestCase):
    def test_count_matches_request(self):
        rng = random.Random(42)
        svs = generate_person_svs(rng, ["22"], 1_000_000, n_svs=7)
        self.assertEqual(len(svs), 7)

    def test_zero_returns_empty(self):
        rng = random.Random(0)
        self.assertEqual(
            generate_person_svs(rng, ["22"], 1_000_000, n_svs=0), [])

    def test_all_records_well_formed(self):
        rng = random.Random(7)
        svs = generate_person_svs(rng, ["22"], 1_000_000, n_svs=20)
        for s in svs:
            self.assertIn(s["svtype"], SV_TYPES)
            self.assertIn(s["alts"][0], {"<DEL>", "<DUP>", "<INV>"})
            self.assertEqual(s["alts"][0], f"<{s['svtype']}>")
            self.assertGreaterEqual(s["pos"], 1)
            self.assertEqual(s["end"], s["pos"] + abs(s["svlen"]))
            self.assertIn(s["gt"], {"0|1", "1|1"})
            if s["svtype"] == "DEL":
                self.assertLess(s["svlen"], 0)
            else:
                self.assertGreater(s["svlen"], 0)

    def test_lengths_within_bounds(self):
        rng = random.Random(1)
        svs = generate_person_svs(rng, ["22"], 1_000_000,
                                  n_svs=50, length_min_bp=100,
                                  length_max_bp=5_000)
        for s in svs:
            self.assertGreaterEqual(abs(s["svlen"]), 100)
            self.assertLessEqual(abs(s["svlen"]), 5_000)

    def test_end_inside_chrom_span(self):
        rng = random.Random(2)
        span = 100_000
        svs = generate_person_svs(rng, ["22"], span,
                                  n_svs=50, length_max_bp=10_000)
        for s in svs:
            self.assertLessEqual(s["end"], span)

    def test_distributes_across_chromosomes(self):
        rng = random.Random(3)
        svs = generate_person_svs(rng, ["19", "20", "21", "22"],
                                  1_000_000, n_svs=200)
        chroms_seen = {s["chrom"] for s in svs}
        # With 200 draws across 4 chroms we should hit every chrom.
        self.assertEqual(chroms_seen, {"19", "20", "21", "22"})

    def test_type_distribution_roughly_matches_weights(self):
        """Default weights are 0.50 / 0.30 / 0.20 — a 1k draw should
        land DEL within ±0.07 of 0.50 and similar for the rest."""
        rng = random.Random(4)
        svs = generate_person_svs(rng, ["22"], 5_000_000, n_svs=1000)
        counts = {t: sum(1 for s in svs if s["svtype"] == t)
                  for t in SV_TYPES}
        total = sum(counts.values())
        del_frac = counts["DEL"] / total
        dup_frac = counts["DUP"] / total
        inv_frac = counts["INV"] / total
        self.assertAlmostEqual(del_frac, 0.50, delta=0.07)
        self.assertAlmostEqual(dup_frac, 0.30, delta=0.07)
        self.assertAlmostEqual(inv_frac, 0.20, delta=0.07)

    def test_reproducible_under_seed(self):
        a = generate_person_svs(random.Random(99), ["22"],
                                1_000_000, n_svs=10)
        b = generate_person_svs(random.Random(99), ["22"],
                                1_000_000, n_svs=10)
        self.assertEqual(len(a), len(b))
        for x, y in zip(a, b):
            self.assertEqual(
                (x["chrom"], x["pos"], x["svtype"], x["svlen"],
                 x["end"], x["gt"], x["ref"]),
                (y["chrom"], y["pos"], y["svtype"], y["svlen"],
                 y["end"], y["gt"], y["ref"]),
            )

    def test_different_seeds_differ(self):
        a = generate_person_svs(random.Random(1), ["22"],
                                1_000_000, n_svs=5)
        b = generate_person_svs(random.Random(2), ["22"],
                                1_000_000, n_svs=5)
        positions_a = [s["pos"] for s in a]
        positions_b = [s["pos"] for s in b]
        self.assertNotEqual(positions_a, positions_b)

    def test_rejects_too_small_chrom(self):
        rng = random.Random(0)
        # Default max length is 10kb — a 5kb chromosome is too small.
        with self.assertRaises(ValueError):
            generate_person_svs(rng, ["22"], 5_000, n_svs=3,
                                length_max_bp=10_000)

    def test_rejects_no_chromosomes(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            generate_person_svs(rng, [], 1_000_000, n_svs=3)

    def test_default_length_constants(self):
        # Lock in the documented default range so a behaviour change is
        # surfaced as a test failure.
        self.assertEqual(DEFAULT_SV_LENGTH_MIN_BP, 50)
        self.assertEqual(DEFAULT_SV_LENGTH_MAX_BP, 10_000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
