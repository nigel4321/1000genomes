"""Tests for the M11 truth-set BED writer."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from syntheticgen.truth import (
    GOLDEN_CATEGORIES,
    TruthBedWriter,
    classify_golden,
    golden_bed_line,
    noise_bed_line,
)


def _snv(**overrides):
    base = {
        "chrom": "chr22",
        "pos": 1000,
        "id": ".",
        "ref": "A",
        "alts": ("G",),
    }
    base.update(overrides)
    return base


class TestClassifyGolden(unittest.TestCase):

    def test_highlighted_wins_over_everything(self):
        v = _snv(id="rs1", clnsig="Pathogenic", cosmic_id="COSV1",
                 svtype="DEL")
        self.assertEqual(classify_golden(v, is_hi=True), "HIGHLIGHTED")

    def test_clinvar_beats_cosmic_and_rsid(self):
        v = _snv(id="rs1", clnsig="Pathogenic", cosmic_id="COSV1")
        self.assertEqual(classify_golden(v, is_hi=False), "CLINVAR")

    def test_cosmic_beats_rsid(self):
        v = _snv(id="rs1", cosmic_id="COSV1")
        self.assertEqual(classify_golden(v, is_hi=False), "COSMIC")

    def test_sv_when_no_other_tags(self):
        v = _snv(svtype="DEL", svlen=-500, end=1500)
        self.assertEqual(classify_golden(v, is_hi=False), "SV")

    def test_rsid_only(self):
        v = _snv(id="rs999")
        self.assertEqual(classify_golden(v, is_hi=False), "RSID")

    def test_dot_clnsig_is_treated_as_missing(self):
        v = _snv(id="rs999", clnsig=".")
        self.assertEqual(classify_golden(v, is_hi=False), "RSID")

    def test_unannotated_record_returns_none(self):
        v = _snv(id=".")
        self.assertIsNone(classify_golden(v, is_hi=False))

    def test_categories_constant_lists_priority_order(self):
        self.assertEqual(GOLDEN_CATEGORIES,
                         ("HIGHLIGHTED", "CLINVAR", "COSMIC", "RSID", "SV"))


class TestGoldenBedLine(unittest.TestCase):

    def test_snv_emits_half_open_interval(self):
        v = _snv(pos=1000, ref="A", alts=("G",), id="rs1")
        line = golden_bed_line(v, "RSID", "0|1")
        chrom, start, end, payload = line.split("\t")
        self.assertEqual(chrom, "chr22")
        self.assertEqual(start, "999")
        self.assertEqual(end, "1000")
        self.assertIn("flag=RSID", payload)
        self.assertIn("id=rs1", payload)
        self.assertIn("gt=0|1", payload)

    def test_deletion_ref_extends_end(self):
        v = _snv(pos=1000, ref="ACGT", alts=("A",))
        line = golden_bed_line(v, "RSID", "0|1")
        _, start, end, _ = line.split("\t")
        self.assertEqual(start, "999")
        self.assertEqual(end, "1003")

    def test_sv_uses_explicit_end(self):
        v = _snv(svtype="DEL", svlen=-500, end=1500)
        line = golden_bed_line(v, "SV", "0|1")
        _, start, end, payload = line.split("\t")
        self.assertEqual(start, "999")
        self.assertEqual(end, "1500")
        self.assertIn("svtype=DEL", payload)
        self.assertIn("svlen=-500", payload)

    def test_clinvar_payload_carries_clnsig(self):
        v = _snv(id="rs1", clnsig="Pathogenic", clndn="Cardiac")
        line = golden_bed_line(v, "CLINVAR", "1|1")
        payload = line.split("\t")[3]
        self.assertIn("flag=CLINVAR", payload)
        self.assertIn("clnsig=Pathogenic", payload)
        self.assertIn("clndn=Cardiac", payload)

    def test_payload_escapes_separators(self):
        v = _snv(id="rs1", clnsig="A;B\tC\nD")
        line = golden_bed_line(v, "CLINVAR", "0|1")
        # Payload field must not contain raw tabs/semicolons in values
        # (semicolons separate keys; tabs separate columns).
        payload = line.split("\t")[3]
        self.assertNotIn("\n", payload)
        # The escaped clnsig should appear with substituted chars
        self.assertIn("clnsig=A,B C D", payload)


class TestNoiseBedLine(unittest.TestCase):

    def test_flip_records_truth_and_called(self):
        v = _snv(id="rs1")
        line = noise_bed_line(v, "FLIP", truth_gt="0|0", called_gt="0|1")
        payload = line.split("\t")[3]
        self.assertIn("flag=FLIP", payload)
        self.assertIn("truth_gt=0|0", payload)
        self.assertIn("called_gt=0|1", payload)

    def test_dropout_records_called_missing(self):
        v = _snv(id="rs1")
        line = noise_bed_line(v, "DROPOUT", "1|1", "./.")
        payload = line.split("\t")[3]
        self.assertIn("flag=DROPOUT", payload)
        self.assertIn("called_gt=./.", payload)


class TestTruthBedWriter(unittest.TestCase):

    def test_writes_two_files_with_sorted_rows(self):
        contig_order = {"chr1": 0, "chr22": 1}
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            golden = tdp / "g.bed"
            noise = tdp / "n.bed"
            w = TruthBedWriter(golden, noise, contig_order=contig_order)
            # Add out-of-order: chr22 first, then chr1
            w.add_golden(_snv(chrom="chr22", pos=2000, id="rs2"),
                         "RSID", "0|1")
            w.add_golden(_snv(chrom="chr1", pos=500, id="rs1"),
                         "RSID", "0|1")
            w.add_golden(_snv(chrom="chr1", pos=100, id="rs0"),
                         "RSID", "0|1")
            w.close()

            lines = golden.read_text().strip().splitlines()
            # contig_order puts chr1 (idx 0) before chr22 (idx 1), and
            # within chr1 the rows sort by pos ascending.
            self.assertEqual(len(lines), 3)
            self.assertTrue(lines[0].startswith("chr1\t99\t"))
            self.assertTrue(lines[1].startswith("chr1\t499\t"))
            self.assertTrue(lines[2].startswith("chr22\t1999\t"))

    def test_counts_track_writes(self):
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            w = TruthBedWriter(tdp / "g.bed", tdp / "n.bed")
            w.add_golden(_snv(id="rs1"), "RSID", "0|1")
            w.add_golden(_snv(id="rs2"), "RSID", "0|1")
            w.add_noise(_snv(id="rs1"), "FLIP", "0|0", "0|1")
            self.assertEqual(w.golden_count, 2)
            self.assertEqual(w.noise_count, 1)
            w.close()

    def test_context_manager_closes(self):
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            golden = tdp / "g.bed"
            noise = tdp / "n.bed"
            with TruthBedWriter(golden, noise) as w:
                w.add_golden(_snv(id="rs1"), "RSID", "0|1")
            self.assertTrue(golden.exists())
            self.assertTrue(noise.exists())

    def test_empty_writer_creates_empty_files(self):
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            golden = tdp / "g.bed"
            noise = tdp / "n.bed"
            TruthBedWriter(golden, noise).close()
            self.assertTrue(golden.exists())
            self.assertEqual(golden.read_text(), "")
            self.assertTrue(noise.exists())
            self.assertEqual(noise.read_text(), "")

    def test_creates_parent_dirs(self):
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            golden = tdp / "nested" / "deeper" / "g.bed"
            noise = tdp / "nested" / "deeper" / "n.bed"
            w = TruthBedWriter(golden, noise)
            w.add_golden(_snv(id="rs1"), "RSID", "0|1")
            w.close()
            self.assertTrue(golden.exists())


if __name__ == "__main__":
    unittest.main()
