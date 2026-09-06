"""Tests for primer_design.py (Step 5)."""
import pytest

from tar_crispr.primer_design import (
    design_tailed_primer,
    design_tailed_primers,
    export_primers,
    _calculate_tm,
    _wallace_tm,
    _gc_percent,
    _check_hairpin,
    _check_dimer,
    _check_gc_clamp,
    _optimize_annealing_length,
    _reverse_complement,
)
from tar_crispr.config import PipelineConfig, HomologyArm


class TestReverseComplement:
    def test_basic(self):
        assert _reverse_complement("ACGT") == "ACGT"
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("ATGC") == "GCAT"


class TestCalculateTm:
    def test_known_tm(self):
        tm = _calculate_tm("ACGTACGTACGTACGTACGT")
        assert 50 < tm < 70  # Reasonable range

    def test_short_seq(self):
        tm = _calculate_tm("ACGTACGT")
        assert tm != 0.0  # Returns a value (may be low for very short seq)

    def test_empty(self):
        assert _calculate_tm("") == 0.0

    def test_wallace(self):
        tm = _wallace_tm("ACGT")
        # Wallace rule for <14bp: 2*(A+T) + 4*(G+C) = 2*2 + 4*2 = 12
        assert tm == 12.0

    def test_wallace_long(self):
        tm = _wallace_tm("ACGTACGTACGTACGTACGT")  # 20bp
        assert tm > 0


class TestGcPercent:
    def test_all_gc(self):
        assert _gc_percent("GCGCGC") == 100.0

    def test_no_gc(self):
        assert _gc_percent("AATATA") == 0.0

    def test_mixed(self):
        assert _gc_percent("GCAT") == 50.0

    def test_empty(self):
        assert _gc_percent("") == 0.0


class TestHairpin:
    def test_no_hairpin(self):
        # A sequence without strong self-complementarity
        result = _check_hairpin("AAAAATTTTTGGGGG")
        # May or may not have hairpin; just testing it runs
        assert isinstance(result, bool)

    def test_with_hairpin(self):
        # Strong palindromic sequence
        seq = "GCGCGCGCGCGCGGCGCGCGCG"
        result = _check_hairpin(seq)
        assert isinstance(result, bool)


class TestDimer:
    def test_no_dimer(self):
        seq = "ACGTACGTACGTACGTACGT"
        other = "TTTTTTTTTTTTTTTTTTTT"
        result = _check_dimer(seq, other)
        assert isinstance(result, bool)

    def test_with_dimer(self):
        seq = "ACGTACGTACGT"
        other = "ACGTACGTACGT"  # Identical -> strong dimer
        result = _check_dimer(seq, other)
        # May detect dimer
        assert isinstance(result, bool)


class TestGcClamp:
    def test_with_clamp(self):
        assert _check_gc_clamp("ACGTACGTACGCGCGC") is True

    def test_without_clamp(self):
        assert _check_gc_clamp("ACGTACGTATAATAAT") is False

    def test_short(self):
        assert _check_gc_clamp("ACGT") is False


class TestOptimizeAnnealingLength:
    def test_finds_in_range(self):
        # Vector sequence with known Tm characteristics
        vec = "GCGCGCGCATATATATATACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        best = _optimize_annealing_length(vec, 10, 30, 18, 25, 58, 62)
        assert 18 <= len(best) <= 25
        assert len(best) > 0

    def test_fallback(self):
        vec = "A" * 100
        best = _optimize_annealing_length(vec, 10, 30, 18, 25, 58, 62)
        assert best is not None
        assert len(best) > 0


class TestDesignTailedPrimer:
    def test_left_primer(self):
        arm = HomologyArm(
            sequence="GCGCGCGCGCATATATATATACGCGCGCGC",
            end="left",
            length=30,
            gc_percent=50.0,
            start_pos=0,
            end_pos=30,
            uniqueness=True,
            secondary_structure=False,
            issues=[],
        )
        # Create a vector with a good annealing region
        vec = "A" * 100 + "GCGCGCGCGCATATATATATGTAC" + "C" * 200
        # Vector cut at position 126 (end of annealing region)
        config = PipelineConfig(blast_available=False)
        primer = design_tailed_primer(arm, vec, 126, "left", config)
        assert primer.name == "FORWARD_left"
        assert primer.tail == arm.sequence
        assert len(primer.sequence) == len(arm.sequence) + len(primer.annealing_region)
        assert isinstance(primer.tm, float)
        assert isinstance(primer.issues, list)

    def test_right_primer(self):
        arm = HomologyArm(
            sequence="GCGCGCGCGCATATATATATACGCGCGCGC",
            end="right",
            length=30,
            gc_percent=50.0,
            start_pos=0,
            end_pos=30,
            uniqueness=True,
            secondary_structure=False,
            issues=[],
        )
        vec = "A" * 100 + "GTACGTATATATGGCGCGCGC" + "C" * 200
        # Vector cut at position 120 (start of annealing region)
        config = PipelineConfig(blast_available=False)
        primer = design_tailed_primer(arm, vec, 120, "right", config)
        assert primer.name == "REVERSE_right"
        assert primer.tail == arm.sequence
        assert isinstance(primer.tm, float)


class TestDesignTailedPrimers:
    def test_both_primers(self):
        left_arm = HomologyArm(
            sequence="GCGCGCGCGCATATATATATAC", end="left", length=22,
            gc_percent=50.0, start_pos=0, end_pos=22,
            uniqueness=True, secondary_structure=False, issues=[],
        )
        right_arm = HomologyArm(
            sequence="TATATATATACGCGCGCGCGCGCGCGCG", end="right", length=28,
            gc_percent=50.0, start_pos=0, end_pos=28,
            uniqueness=True, secondary_structure=False, issues=[],
        )
        vec = "A" * 50 + "GTACGTACGTACGTACGTACGTACGTACGTACGTACGT" + "C" * 100
        config = PipelineConfig(blast_available=False)
        primers = design_tailed_primers(
            left_arm, right_arm, vec, 50, 90, config
        )
        assert "left" in primers
        assert "right" in primers
        assert len(primers["left"].sequence) > 0
        assert len(primers["right"].sequence) > 0


class TestExportPrimers:
    def test_csv_export(self, tmp_path):
        left_arm = HomologyArm(
            sequence="GCGCGCGCGCATATATATATAC", end="left", length=22,
            gc_percent=50.0, start_pos=0, end_pos=22,
            uniqueness=True, secondary_structure=False, issues=[],
        )
        right_arm = HomologyArm(
            sequence="TATATATATACGCGCGCGCGCG", end="right", length=22,
            gc_percent=50.0, start_pos=0, end_pos=22,
            uniqueness=True, secondary_structure=False, issues=[],
        )
        vec = "A" * 50 + "GTACGTACGTACGTACGTACGTACGTACGTACGTACGT" + "C" * 100
        config = PipelineConfig(blast_available=False)
        primers = design_tailed_primers(left_arm, right_arm, vec, 50, 90, config)
        csv_path = str(tmp_path / "primers.csv")
        result_path = export_primers(primers, csv_path)
        assert result_path == csv_path
        with open(csv_path) as fh:
            content = fh.read()
        assert "Primer_Name" in content
        assert "FORWARD_left" in content
        assert "REVERSE_right" in content
