"""Tests for homology_arms.py (Step 4)."""
import pytest

from tar_crispr.homology_arms import (
    extract_homology_arm,
    validate_gc,
    check_uniqueness,
    check_secondary_structure,
    check_restriction_sites,
    validate_arm,
    find_valid_arm,
    design_homology_arms,
)
from tar_crispr.config import PipelineConfig


class TestExtractHomologyArm:
    def test_left_arm(self):
        frag = "A" * 100
        arm = extract_homology_arm(frag, "left", 50)
        assert arm == "A" * 50

    def test_right_arm(self):
        frag = "A" * 100
        arm = extract_homology_arm(frag, "right", 50)
        assert arm == "A" * 50

    def test_short_fragment(self):
        frag = "ACGT"
        arm = extract_homology_arm(frag, "left", 50)
        assert arm == "ACGT"  # Returns whole fragment if shorter

    def test_invalid_end(self):
        with pytest.raises(ValueError):
            extract_homology_arm("AAAA", "middle", 50)


class TestValidateGC:
    def test_good_gc(self):
        config = PipelineConfig(min_gc=0.40, max_gc=0.65, max_gc_warning=0.75)
        passed, issues = validate_gc("GCGCGCGCATATATATAT", config)
        # GC% = 9/17 = 52.9%
        assert passed
        assert issues == []

    def test_low_gc(self):
        config = PipelineConfig(min_gc=0.40, max_gc=0.65, max_gc_warning=0.75)
        passed, issues = validate_gc("ATATATATATAT", config)
        assert passed is False
        assert any("below minimum" in i for i in issues)

    def test_high_gc(self):
        config = PipelineConfig(min_gc=0.40, max_gc=0.65, max_gc_warning=0.75)
        passed, issues = validate_gc("GCGCGCGCGCGCGCGCGC", config)
        assert passed is False
        assert any("exceeds" in i or "above" in i for i in issues)

    def test_very_high_gc_warning(self):
        config = PipelineConfig(min_gc=0.40, max_gc=0.65, max_gc_warning=0.75)
        passed, issues = validate_gc("G" * 20, config)
        assert passed is False
        assert any("exceeds warning" in i for i in issues)


class TestCheckUniqueness:
    def test_unique_in_genome(self):
        config = PipelineConfig()
        arm = "ACGTACGTACGTACGTACGT"
        genome = "TTTT" + arm + "GGGG"
        vector = "CCCC" * 100
        passed, issues = check_uniqueness(arm, genome, vector, config)
        assert passed
        assert issues == []

    def test_non_unique(self):
        config = PipelineConfig()
        arm = "ACGTACGTACGTACGTACGT"
        genome = arm + "NNNN" + arm  # appears twice
        vector = "CCCC" * 100
        passed, issues = check_uniqueness(arm, genome, vector, config)
        assert passed is False
        assert any("non-unique" in i for i in issues)

    def test_in_vector(self):
        config = PipelineConfig()
        arm = "ACGTACGTACGTACGTACGT"
        genome = "TTTT" + arm + "GGGG"
        vector = arm + "NNNN" * 10
        passed, issues = check_uniqueness(arm, genome, vector, config)
        assert passed is False
        assert any("capture vector" in i for i in issues)


class TestCheckSecondaryStructure:
    def test_no_structure(self):
        passed, issues = check_secondary_structure("AAAAATTTTTGGGGGCCCCC")
        assert passed
        assert issues == []

    def test_with_inverted_repeat(self):
        # Create a sequence with a 6bp inverted repeat
        seq = "ACGTAC" + "TTTTTT" + "GTACGT"
        # "ACGTAC" at start, "GTACGT" at end -> RC of "ACGTAC" is "GTACGT"
        passed, issues = check_secondary_structure(seq)
        assert passed is False
        assert any("inverted" in i.lower() for i in issues)


class TestCheckRestrictionSites:
    def test_no_enzymes(self):
        passed, issues = check_restriction_sites("ACGT", [])
        assert passed

    def test_site_present(self):
        passed, issues = check_restriction_sites("GAATTC", ["EcoRI"])
        assert passed is False
        assert any("EcoRI" in i for i in issues)

    def test_site_absent(self):
        passed, issues = check_restriction_sites("AAAA", ["EcoRI"])
        assert passed
        assert issues == []

    def test_unknown_enzyme(self):
        passed, issues = check_restriction_sites("AAAA", ["FakeEnz"])
        assert passed is False
        assert any("Unknown" in i for i in issues)


class TestValidateArm:
    def test_valid_arm(self):
        config = PipelineConfig(min_gc=0.40, max_gc=0.65, max_gc_warning=0.75)
        arm = "GCGCGCGCATATATATATAC"  # GC% = 55%
        genome = "TTTT" + arm + "GGGG" * 20
        vector = "CCCC" * 100
        passed, issues = validate_arm(arm, genome, vector, config)
        assert passed


class TestFindValidArm:
    def test_finds_valid_arm(self):
        config = PipelineConfig(
            homology_arm_length=20,
            min_gc=0.40, max_gc=0.65, max_gc_warning=0.75,
            max_shift=20, shift_increment=5,
        )
        # Fragment: 20bp good arm + 20bp bad arm + 100bp middle
        good_arm = "GCGCGCGCATATATATATAC"
        bad_arm = "TTTTTTTTTTTTTTTTTTTT"
        frag = good_arm + bad_arm + "A" * 100
        genome = "AAAA" + frag + "TTTT"
        vector = "CCCC" * 100
        arm = find_valid_arm(frag, "left", config, genome, vector)
        assert arm is not None
        assert arm.length == 20
        assert arm.sequence == good_arm

    def test_shifts_for_valid_arm(self):
        config = PipelineConfig(
            homology_arm_length=20,
            min_gc=0.40, max_gc=0.65, max_gc_warning=0.75,
            max_shift=40, shift_increment=5,
        )
        # First 25bp are AT-rich (GC=0%), first valid 20bp window with GC>=40%
        # only appears after shift=25
        bad_first = "AT" * 13  # 26bp of AT (GC=0%)
        good_arm = "GCGCGCGCATATATATATAC"  # GC=55%
        frag = bad_first + good_arm + "A" * 100
        genome = "GGGG" + frag + "CCCC"
        vector = "TTTT" * 100
        arm = find_valid_arm(frag, "left", config, genome, vector)
        assert arm is not None
        # Should have found the good_arm (first window with GC >= 40%)
        # At shift=10: GC = 6/20 = 30% (fails)
        # At shift=15: GC = 6/20 = 30% (fails)
        # At shift=20: "AT" + 18bp of good_arm, GC = 10/20 = 50% (passes!)
        # At shift=25: good_arm, GC = 11/20 = 55% (passes)
        # The algorithm finds first valid one, which is shift=20
        assert arm.gc_percent >= 40.0
        assert arm.length == 20

    def test_best_when_none_pass(self):
        config = PipelineConfig(
            homology_arm_length=20,
            min_gc=0.40, max_gc=0.65, max_gc_warning=0.75,
            max_shift=10, shift_increment=5,
        )
        # All arms will have issues — should return best attempt
        frag = "GC" * 60  # Very high GC, will fail
        genome = "AAAA" + frag + "TTTT"
        vector = "CCCC" * 100
        arm = find_valid_arm(frag, "left", config, genome, vector)
        assert arm is not None  # Returns best attempt


class TestDesignHomologyArms:
    def test_both_arms(self):
        config = PipelineConfig(
            homology_arm_length=20,
            min_gc=0.0,  # Relax GC for testing
            max_gc=1.0,
            max_gc_warning=1.0,
            max_shift=20,
            shift_increment=5,
        )
        arm_seq = "GCGCGCGCGCATATATATAT"
        frag = arm_seq + "A" * 60 + reverse_comp(arm_seq) if False else arm_seq + "A" * 60 + arm_seq
        genome = "AAAA" + frag + "TTTT"
        vector = "CCCC" * 100
        result = design_homology_arms(frag, genome, vector, config)
        assert "left" in result
        assert "right" in result
        assert result["left"].sequence is not None
        assert result["right"].sequence is not None


def reverse_comp(seq):
    from tar_crispr.pam_finder import reverse_complement
    return reverse_complement(seq)
