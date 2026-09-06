"""Tests for fragment_ends.py (Step 3)."""
import pytest

from tar_crispr.fragment_ends import (
    FragmentEnds,
    compute_cut_position,
    extract_fragment,
    get_cut_junction_sequences,
)
from tar_crispr.config import PAMCandidate


class TestComputeCutPosition:
    def test_forward_cut(self):
        sg = PAMCandidate(position=100, strand="+", protospacer="A" * 20,
                          pam="NGG", cut_position=117, gc_percent=0.0,
                          polyt_flag=False)
        assert compute_cut_position(sg, "left") == 117

    def test_reverse_cut(self):
        sg = PAMCandidate(position=500, strand="-", protospacer="A" * 20,
                          pam="NGG", cut_position=523, gc_percent=0.0,
                          polyt_flag=False)
        assert compute_cut_position(sg, "right") == 523


class TestExtractFragment:
    def test_basic_extraction(self):
        full = "A" * 100 + "B" * 200 + "C" * 100
        left_sg = PAMCandidate(position=0, strand="+", protospacer="A" * 20,
                               pam="NGG", cut_position=100, gc_percent=0.0,
                               polyt_flag=False)
        right_sg = PAMCandidate(position=0, strand="+", protospacer="C" * 20,
                                pam="NGG", cut_position=300, gc_percent=0.0,
                                polyt_flag=False)
        frag = extract_fragment(full, left_sg, right_sg)
        assert frag.fragment_length == 200
        assert frag.sequence == "B" * 200
        assert frag.left_cut_pos == 100
        assert frag.right_cut_pos == 300

    def test_left_after_right_raises(self):
        full = "A" * 1000
        left_sg = PAMCandidate(position=0, strand="+", protospacer="A" * 20,
                               pam="NGG", cut_position=500, gc_percent=0.0,
                               polyt_flag=False)
        right_sg = PAMCandidate(position=0, strand="+", protospacer="C" * 20,
                                pam="NGG", cut_position=100, gc_percent=0.0,
                                polyt_flag=False)
        with pytest.raises(ValueError, match="before right cut"):
            extract_fragment(full, left_sg, right_sg)

    def test_gc_content(self):
        full = "A" * 100 + "GCGCGC" + "C" * 100
        left_sg = PAMCandidate(position=0, strand="+", protospacer="A" * 20,
                               pam="NGG", cut_position=100, gc_percent=0.0,
                               polyt_flag=False)
        right_sg = PAMCandidate(position=0, strand="+", protospacer="C" * 20,
                                pam="NGG", cut_position=106, gc_percent=0.0,
                                polyt_flag=False)
        frag = extract_fragment(full, left_sg, right_sg)
        assert frag.gc_percent == 100.0  # All GC

    def test_cut_junction_sequences(self):
        full = "A" * 1000
        left_sg = PAMCandidate(position=0, strand="+", protospacer="A" * 20,
                               pam="NGG", cut_position=100, gc_percent=0.0,
                               polyt_flag=False)
        upstream, downstream = get_cut_junction_sequences(left_sg, full, 20)
        assert upstream == "A" * 20
        assert downstream == "A" * 20
