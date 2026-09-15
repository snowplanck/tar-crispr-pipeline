"""Regression tests for bugs found and fixed during real-data validation
(BGC de colinomicina, vetor pCAP03). Each test pins a specific failure
mode so it cannot silently reappear.
"""
from tar_crispr.config import PipelineConfig, HomologyArm
from tar_crispr.primer_design import (
    design_tailed_primer,
    design_tailed_primers,
    _optimize_annealing_length,
)


# ---------------------------------------------------------------------------
# Bug: _optimize_annealing_length used to slide past the cut site (positive
# shift), landing the annealing region on the WRONG side of the cut.
# ---------------------------------------------------------------------------
def test_optimize_annealing_never_crosses_cut_downstream():
    # AT-rich desert right at the cut, GC-rich region far upstream.
    at_desert = "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA"
    gc_rich_upstream = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
    vector = gc_rich_upstream + at_desert
    cut = len(vector) - 5  # cut sits inside the AT desert

    result = _optimize_annealing_length(
        vector, cut,
        min_len=18, max_len=25,
        tm_min=58.0, tm_max=62.0,
        max_shift=100,
    )

    found_at = vector.find(result)
    assert found_at != -1, "annealing region must be a real substring of the vector"
    end_of_region = found_at + len(result)
    # The 3' end of the annealing region must not be placed AFTER the cut
    # (i.e. must not have slid downstream past it).
    assert end_of_region <= cut + 3, (
        f"annealing region ends at {end_of_region}, more than 3bp past the "
        f"cut at {cut} — regression of the 'wrong side of cut' bug"
    )


# ---------------------------------------------------------------------------
# Bug: self-dimer (primer vs itself) and cross-dimer (forward vs reverse)
# were conflated into a single 'dimer' field, so a primer with a
# self-complementary 3' end got flagged as a fatal cross-dimer issue.
# ---------------------------------------------------------------------------
def test_self_dimer_and_cross_dimer_are_independent_fields():
    config = PipelineConfig()
    vector = "ACGT" * 40 + "GTTTAAAC" + "ACGT" * 40  # SwaI-like site in the middle
    cut = vector.find("GTTTAAAC") + 4

    left_arm = HomologyArm(
        sequence="C" * 30 + "GCGC",  # palindromic 3' end -> self-dimer prone
        end="left", length=34, gc_percent=70.0,
        start_pos=0, end_pos=34, uniqueness=True, secondary_structure=False,
    )
    right_arm = HomologyArm(
        sequence="T" * 30 + "AAAA",
        end="right", length=34, gc_percent=10.0,
        start_pos=0, end_pos=34, uniqueness=True, secondary_structure=False,
    )

    primers = design_tailed_primers(left_arm, right_arm, vector, cut, cut, config)
    left, right = primers["left"], primers["right"]

    # Both fields must exist independently (regression: 'dimer' attribute
    # no longer exists — this raises AttributeError if the split broke).
    assert hasattr(left, "self_dimer")
    assert hasattr(left, "cross_dimer")
    assert hasattr(right, "self_dimer")
    assert hasattr(right, "cross_dimer")
    assert isinstance(left.cross_dimer, bool)
    assert isinstance(right.cross_dimer, bool)
