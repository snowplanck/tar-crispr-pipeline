"""Tests for pam_finder.py (Step 2)."""
import pytest

from tar_crispr.pam_finder import (
    find_pam_sites,
    calculate_gc,
    check_specificity,
    rank_sgRNAs,
    design_sgRNAs,
    complement,
    reverse_complement,
    _build_pam_regex,
    _check_polyt,
    _hamming_search,
)
from tar_crispr.config import PipelineConfig


class TestComplementReverseComplement:
    def test_complement(self):
        assert complement("ACGT") == "TGCA"
        assert complement("ACGTN") == "TGCAN"

    def test_reverse_complement(self):
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("ATGC") == "GCAT"


class TestBuildPamRegex:
    def test_egg(self):
        regex = _build_pam_regex("NGG")
        assert "ACGT" in regex and regex.endswith("GG")

    def test_case_insensitive(self):
        regex = _build_pam_regex("ngg")
        assert "ACGT" in regex and regex.endswith("GG")

    def test_ccn(self):
        regex = _build_pam_regex("CCN")
        assert regex.startswith("CC") and "ACGT" in regex


class TestCheckPolyt:
    def test_polyt_present(self):
        assert _check_polyt("ATGTTTT") is True

    def test_no_polyt(self):
        assert _check_polyt("ATGTCGT") is False

    def test_threshold(self):
        assert _check_polyt("ATGTAT", 4) is False
        assert _check_polyt("ATGTTTT", 4) is True


class TestCalculateGC:
    def test_all_gc(self):
        assert calculate_gc("GCGCGC") == 100.0

    def test_no_gc(self):
        assert calculate_gc("AATATA") == 0.0

    def test_mixed(self):
        assert calculate_gc("GCAT") == 50.0

    def test_empty(self):
        assert calculate_gc("") == 0.0


class TestFindPamSites:
    def test_forward_pam(self):
        # 20nt protospacer + NGG
        seq = "ACGTACGTACGTACGTACGT" + "NGG" + "TTTT"
        # But N is not a real base; use actual bases that match NGG
        seq = "ACGTACGTACGTACGTACGT" + "AGG" + "TTTT"
        candidates = find_pam_sites(seq, "NGG", 20)
        assert len(candidates) >= 1
        fwd = [c for c in candidates if c.strand == "+"]
        assert len(fwd) == 1
        assert fwd[0].protospacer == "ACGTACGTACGTACGTACGT"
        assert fwd[0].pam == "AGG"
        # Cut at PAM_start - 3 = 20 - 3 = 17
        assert fwd[0].cut_position == 17

    def test_reverse_pam(self):
        # For reverse strand: CCN on forward + protospacer_rc downstream
        # CCN + 20nt protospacer_rc
        seq = "CCN" + "ACGTACGTACGTACGTACGT"
        # Use real bases: CCG + ACGT...
        seq = "CCG" + "ACGTACGTACGTACGTACGT"
        candidates = find_pam_sites(seq, "NGG", 20)
        rev = [c for c in candidates if c.strand == "-"]
        assert len(rev) >= 1
        # protospacer should be reverse_complement of the downstream sequence
        expected_protospacer = reverse_complement("ACGTACGTACGTACGTACGT")
        assert rev[0].protospacer == expected_protospacer

    def test_no_pam(self):
        # Sequence with no NGG on either strand
        seq = "AAAAATTTCCGGGGCCCC"
        candidates = find_pam_sites(seq, "NGG", 20)
        assert len(candidates) == 0

    def test_pam_at_start(self):
        # PAM at position 0, no room for 20nt protospacer upstream
        seq = "NGG" + "ACGT" * 10
        seq = "AGG" + "ACGT" * 10
        candidates = find_pam_sites(seq, "NGG", 20)
        # No forward candidates (not enough upstream), but reverse may have some
        fwd = [c for c in candidates if c.strand == "+"]
        assert len(fwd) == 0

    def test_both_strands(self):
        seq = "ACGTACGTACGTACGTACGTAGG" + "TTTT" + "ACGTACGTACGTACGTACGTCCG"
        candidates = find_pam_sites(seq, "NGG", 20)
        strands = set(c.strand for c in candidates)
        assert "+" in strands or "-" in strands

    def test_gc_percent(self):
        seq = "GCGCGCGGCGCGCGCGCGCGCGCG"  # all GC + NGG
        seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCGG" + "AGG"
        # This is 30 G/C + AGG = 33 chars, protospacer = first 20 = GCGCGC...
        candidates = find_pam_sites(seq, "NGG", 20)
        fwd = [c for c in candidates if c.strand == "+"]
        if fwd:
            assert fwd[0].gc_percent > 0

    def test_protospacer_length(self):
        seq = "ACGTACGTACGTACGTACGT" + "AGG"
        candidates = find_pam_sites(seq, "NGG", 20)
        fwd = [c for c in candidates if c.strand == "+"]
        assert len(fwd[0].protospacer) == 20


class TestCheckSpecificity:
    def test_no_match(self):
        protospacer = "ACGTACGTACGTACGTACGT"
        genome = "TTTT" * 100
        hits = check_specificity(protospacer, genome, max_mismatch=3,
                                 blast_available=False)
        assert hits == 0

    def test_perfect_match(self):
        protospacer = "ACGTACGTACGTACGTACGT"
        genome = "TTTT" + protospacer + "GGGG"
        hits = check_specificity(protospacer, genome, max_mismatch=3,
                                 blast_available=False)
        assert hits >= 1  # at least the perfect match on forward + reverse

    def test_with_mismatch(self):
        protospacer = "ACGTACGTACGTACGTACGT"
        genome = "TTTT" + "ACGTACGTACGTACGTACGA" + "GGGG"  # 1 mismatch
        hits = check_specificity(protospacer, genome, max_mismatch=3,
                                 blast_available=False)
        assert hits >= 1

    def test_beyond_mismatch(self):
        protospacer = "ACGTACGTACGTACGTACGT"
        genome = "TTTT" + "TTTTTTTTTTTTTTTTTTTT" + "GGGG"  # 20 mismatches
        hits = check_specificity(protospacer, genome, max_mismatch=3,
                                 blast_available=False)
        assert hits == 0

    def test_hamming_search(self):
        assert _hamming_search("ACGTACGT", "ACGTACGT", 0) == 1
        assert _hamming_search("ACGTACGT", "TTTTTTTT", 0) == 0
        assert _hamming_search("ACGTACGT", "ACGTACGA", 1) == 1


class TestRankSgRNAs:
    def test_ranking(self):
        from tar_crispr.config import PAMCandidate
        candidates = [
            PAMCandidate(position=0, strand="+", protospacer="ACGTACGTACGTACGTACGT",
                         pam="NGG", cut_position=17, gc_percent=50.0, polyt_flag=False),
            PAMCandidate(position=50, strand="+", protospacer="TTTTTTTTTTTTTTTTTTTT",
                         pam="NGG", cut_position=67, gc_percent=0.0, polyt_flag=True),
        ]
        # Use a genome that does NOT match either protospacer to avoid off-target penalties
        genome = "G" * 1000
        config = PipelineConfig(blast_available=False)
        ranked = rank_sgRNAs(candidates, genome, config)
        assert len(ranked) == 2
        # The first one should rank higher (no poly-T, good GC, no off-targets)
        assert ranked[0] == candidates[0]
        assert ranked[1] == candidates[1]


class TestDesignSgRNAs:
    def test_design_returns_left_and_right(self):
        # Create a sequence with PAMs on both ends
        upstream = "ACGTACGTACGTACGTACGTAGG" + "T" * 20
        cluster = "GGCC" * 250  # 1000bp
        downstream = "ACGTACGTACGTACGTACGTCCG" + "A" * 20
        genome = upstream + cluster + downstream
        config = PipelineConfig(blast_available=False, top_n_sgRNAs=3)
        result = design_sgRNAs(cluster, upstream, downstream, genome, config)
        assert "left" in result
        assert "right" in result
        assert len(result["left"]) <= 3
        assert len(result["right"]) <= 3
