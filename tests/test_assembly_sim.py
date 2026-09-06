"""Tests for assembly_sim.py (Step 6)."""
import pytest

from tar_crispr.assembly_sim import (
    AssemblyDesign,
    simulate_assembly,
    simulate_pydna_assembly,
    _check_junctions,
    _check_spurious_overlaps,
    _count_overlapping,
)
from tar_crispr.config import HomologyArm, AssemblyResult


def make_arm(seq, end="left"):
    return HomologyArm(
        sequence=seq, end=end, length=len(seq),
        gc_percent=50.0, start_pos=0, end_pos=len(seq),
        uniqueness=True, secondary_structure=False, issues=[],
    )


class TestCountOverlapping:
    def test_single(self):
        assert _count_overlapping("AAAAACCCCC", "AAAAA") == 1

    def test_multiple(self):
        assert _count_overlapping("AAAA", "AA") == 3

    def test_overlapping(self):
        assert _count_overlapping("GGGGGG", "GGGG") == 3


class TestCheckJunctions:
    def test_matching_junctions(self):
        issues = []
        fragment = "AAAACGTACG" + "MIDDLE" + "TTT"
        left_arm = make_arm("AAAACGTACG", "left")
        left_arm.end_pos = len(left_arm.sequence)  # start=0, end=10
        right_arm = make_arm("TTT", "right")
        right_arm.start_pos = len(fragment) - 3
        right_arm.end_pos = len(fragment)
        ok = _check_junctions(fragment, left_arm, right_arm, issues)
        assert ok
        assert issues == []

    def test_mismatched_left(self):
        issues = []
        fragment = "TTTTCGTACG" + "MIDDLE" + "TTT"
        left_arm = make_arm("AAAACGTACG", "left")
        right_arm = make_arm("TTT", "right")
        right_arm.start_pos = len(fragment) - 3
        right_arm.end_pos = len(fragment)
        ok = _check_junctions(fragment, left_arm, right_arm, issues)
        assert not ok
        assert any("Left" in i or "left" in i for i in issues)

    def test_short_fragment(self):
        """Fragment too short to contain both arms."""
        issues = []
        left_arm = make_arm("A" * 20, "left")
        right_arm = make_arm("T" * 20, "right")
        fragment = "SHORT"  # 5bp, too short for 40bp of arms
        right_arm.start_pos = 3
        right_arm.end_pos = 5
        ok = _check_junctions(fragment, left_arm, right_arm, issues)
        assert not ok
        assert any("short" in i.lower() for i in issues)

    def test_positions_beyond_fragment(self):
        issues = []
        left_arm = make_arm("A" * 20, "left")
        right_arm = make_arm("T" * 20, "right")
        right_arm.start_pos = 1000  # Beyond fragment
        right_arm.end_pos = 1020
        ok = _check_junctions("ACGT", left_arm, right_arm, issues)
        assert not ok
        assert any("beyond" in i for i in issues)


class TestCheckSpuriousOverlaps:
    def test_no_spurious(self):
        issues = _check_spurious_overlaps("XYZABCDWVUP", "ABCD", "WVUP")
        assert issues == []

    def test_spurious_left(self):
        arm = "AAAA"
        seq = arm + "MID" + arm + arm + "END"
        issues = _check_spurious_overlaps(seq, arm, "TTTT")
        assert len(issues) > 0
        assert any("Left" in i for i in issues)


class TestSimulateAssembly:
    def test_successful_assembly(self):
        """Assembly with matching arms should succeed."""
        left_arm_seq = "GCGCGCGCGCATATATAT"  # 18bp
        right_arm_seq = "TATATATATACGCGCGCG"
        fragment = left_arm_seq + "BGC_SEQUENCE_HERE_100bp" + right_arm_seq
        vector = "VEC5" * 50 + left_arm_seq + "MIDDLE" + right_arm_seq + "VEC3" * 50
        # Actually, for the simulation, the backbone should NOT contain the arms
        # The vector backbone is what remains after cutting out the insert site
        # In TAR, the vector is linearized and the fragment recombines
        # Let's set up: vector has a cut site, backbone = vec[:cut_left] + vec[cut_right:]
        # The fragment has arms that match the backbone ends

        # Simpler approach: vector backbone ends match the arms
        vec_5p = "BACKBONE_LEFT" + left_arm_seq
        vec_3p = right_arm_seq + "BACKBONE_RIGHT"
        vector = vec_5p + "TO_BE_CUT_OUT" + vec_3p
        cut_left = len(vec_5p) - len(left_arm_seq)  # Just before the arm
        cut_right = len(vec_5p) + len("TO_BE_CUT_OUT")  # After the cut-out region

        left_arm = make_arm(left_arm_seq, "left")
        right_arm = make_arm(right_arm_seq, "right")

        design = AssemblyDesign(
            fragment_seq=fragment,
            vector_seq=vector,
            vector_cut_left=cut_left,
            vector_cut_right=cut_right,
            left_arm=left_arm,
            right_arm=right_arm,
        )

        result = simulate_assembly(design)
        assert isinstance(result, AssemblyResult)
        assert result.circular is True
        assert result.final_size > 0

    def test_bad_cut_positions(self):
        left_arm = make_arm("AAAA", "left")
        right_arm = make_arm("TTTT", "right")
        design = AssemblyDesign(
            fragment_seq="AAAAMIDTTTT",
            vector_seq="VEC",
            vector_cut_left=100,
            vector_cut_right=50,  # left > right
            left_arm=left_arm,
            right_arm=right_arm,
        )
        result = simulate_assembly(design)
        assert result.success is False
        assert any("before" in i for i in result.issues)

    def test_pydna_wrapper(self):
        left_arm_seq = "GCGCGCGCGCATATATAT"
        right_arm_seq = "TATATATATACGCGCGCG"
        fragment = left_arm_seq + "MIDDLE_SEQ" + right_arm_seq
        vector = "BACKBONE_LEFT" + left_arm_seq + "CUT_REGION" + right_arm_seq + "BACKBONE_RIGHT"
        cut_left = len("BACKBONE_LEFT" + left_arm_seq) - len(left_arm_seq)
        cut_right = len("BACKBONE_LEFT") + len(left_arm_seq) + len("CUT_REGION")

        left_arm = make_arm(left_arm_seq, "left")
        result = simulate_pydna_assembly(
            fragment, vector, left_arm, right_arm_seq, cut_left, cut_right
        )
        assert isinstance(result, AssemblyResult)

    def test_expected_size(self):
        left_arm_seq = "GCGCGCGCATATATATAT"
        right_arm_seq = "TATATATATACGCGCGCG"
        middle = "BGC"
        fragment = left_arm_seq + middle + right_arm_seq
        vec_5p = "LEFT_BACKBONE" + left_arm_seq
        vec_3p = right_arm_seq + "RIGHT_BACKBONE"
        vector = vec_5p + "CUT_REGION" + vec_3p
        cut_left = len(vec_5p)
        cut_right = len(vec_5p) + len("CUT_REGION")

        # Set proper start/end positions for arms within the fragment
        left_arm = HomologyArm(
            sequence=left_arm_seq, end="left", length=len(left_arm_seq),
            gc_percent=50.0, start_pos=0, end_pos=len(left_arm_seq),
            uniqueness=True, secondary_structure=False, issues=[],
        )
        right_start = len(fragment) - len(right_arm_seq)
        right_arm = HomologyArm(
            sequence=right_arm_seq, end="right", length=len(right_arm_seq),
            gc_percent=50.0, start_pos=right_start, end_pos=len(fragment),
            uniqueness=True, secondary_structure=False, issues=[],
        )

        expected = len(vec_5p) + len(fragment) + len(vec_3p)
        design = AssemblyDesign(
            fragment_seq=fragment,
            vector_seq=vector,
            vector_cut_left=cut_left,
            vector_cut_right=cut_right,
            left_arm=left_arm,
            right_arm=right_arm,
            expected_size=expected,
        )

        result = simulate_assembly(design)
        assert result.final_size == expected
