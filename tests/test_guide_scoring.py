"""Tests for position-weighted, genome-wide sgRNA specificity scoring."""
import random

import pytest

from tar_crispr.config import PAMCandidate, PipelineConfig
from tar_crispr.guide_scoring import (
    MIT_WEIGHTS, PAM_FACTORS, annotate_specificity, mismatch_positions,
    offtarget_hit_score, scan_offtargets, specificity_from_hits,
)
from tar_crispr.pam_finder import (
    design_sgRNAs, find_pam_sites, rank_sgRNAs, reverse_complement,
)

GUIDE = "GACGTTGCAAGCTTGCAGTC"      # 20 nt


def rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def mutate(seq, positions):
    out = list(seq)
    for p in positions:
        out[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[out[p]]
    return "".join(out)


def genome_with(site, pam, strand="+", flank=300, seed=1):
    """Random genome with *site*+*pam* planted (optionally on the reverse strand)."""
    unit = site + pam
    if strand == "-":
        unit = reverse_complement(unit)
    return rnd(flank, seed) + unit + rnd(flank, seed + 100), flank


class TestHitScore:
    def test_perfect_copy_is_one(self):
        assert offtarget_hit_score([]) == 1.0

    def test_single_seed_mismatch_much_lower_than_distal(self):
        assert offtarget_hit_score([19]) < 0.5
        assert offtarget_hit_score([0]) == 1.0
        assert offtarget_hit_score([19]) < offtarget_hit_score([2]) < offtarget_hit_score([0]) + 1e-9

    def test_two_mismatch_formula(self):
        # positions 5 and 15: (1-.395)(1-.828) * 1/((19-10)/19*4+1) * 1/4
        expected = (1 - 0.395) * (1 - 0.828) / (((19 - 10) / 19) * 4 + 1) / 4
        assert offtarget_hit_score([5, 15]) == pytest.approx(expected)

    def test_clustered_mismatches_score_lower_than_spread(self):
        assert offtarget_hit_score([10, 11]) < offtarget_hit_score([2, 18]) or \
            offtarget_hit_score([10, 11]) < offtarget_hit_score([5, 17])

    def test_more_mismatches_lower_score(self):
        assert offtarget_hit_score([5, 8, 12]) < offtarget_hit_score([5, 8])

    def test_weights_length(self):
        assert len(MIT_WEIGHTS) == 20

    def test_mismatch_positions(self):
        assert mismatch_positions("AAAA", "ATAT") == (1, 3)


class TestSpecificityCombination:
    def test_no_hits_is_100(self):
        assert specificity_from_hits([]) == 100.0

    def test_one_perfect_copy_is_50(self):
        assert specificity_from_hits([1.0]) == pytest.approx(50.0)

    def test_monotonic(self):
        assert specificity_from_hits([0.2, 0.1]) < specificity_from_hits([0.2])


class TestScan:
    def _find(self, guide, genome, **kw):
        return scan_offtargets([guide], genome, 3, 12, **kw)[guide]

    def test_finds_perfect_site_forward(self):
        g, flank = genome_with(GUIDE, "AGG")
        hits = self._find(GUIDE, g)
        assert len(hits) == 1
        assert hits[0].cut_position == flank + 17 and hits[0].strand == "+"
        assert hits[0].mismatches == 0 and hits[0].pam_class == "NGG"

    def test_finds_perfect_site_reverse_and_matches_pam_finder(self):
        g, flank = genome_with(GUIDE, "TGG", strand="-")
        hits = self._find(GUIDE, g)
        assert len(hits) == 1 and hits[0].strand == "-"
        finder = [c for c in find_pam_sites(g) if c.strand == "-" and c.protospacer == GUIDE]
        assert hits[0].cut_position == finder[0].cut_position

    def test_distal_mismatches_found_with_positions(self):
        site = mutate(GUIDE, [1, 4])
        g, _ = genome_with(site, "CGG")
        h = self._find(GUIDE, g)[0]
        assert h.mismatch_positions == (1, 4) and h.mismatches == 2

    def test_one_seed_mismatch_found_two_seed_mismatches_not(self):
        one = mutate(GUIDE, [15])
        two = mutate(GUIDE, [15, 18])
        assert self._find(GUIDE, genome_with(one, "AGG")[0])
        assert not self._find(GUIDE, genome_with(two, "AGG")[0])

    def test_over_tolerance_ignored(self):
        site = mutate(GUIDE, [0, 1, 3, 4])
        assert not self._find(GUIDE, genome_with(site, "AGG")[0])

    def test_pam_required(self):
        g, _ = genome_with(GUIDE, "ATT")
        assert self._find(GUIDE, g) == []

    def test_noncanonical_pam_class(self):
        assert self._find(GUIDE, genome_with(GUIDE, "AAG")[0])[0].pam_class == "NAG"
        assert self._find(GUIDE, genome_with(GUIDE, "AGA")[0])[0].pam_class == "NGA"

    def test_multiple_guides_one_pass(self):
        g2 = "GTTGCAAGCTTGCAGTCACG"
        genome = genome_with(GUIDE, "AGG", seed=1)[0] + genome_with(g2, "TGG", seed=2)[0]
        res = scan_offtargets([GUIDE, g2], genome, 3, 12)
        assert len(res[GUIDE]) == 1 and len(res[g2]) == 1


class TestAnnotate:
    def _cand(self, genome, guide, flank, strand="+"):
        c = [c for c in find_pam_sites(genome) if c.protospacer == guide and c.strand == strand]
        assert c
        return c[0]

    def test_unique_guide_scores_100(self):
        g, flank = genome_with(GUIDE, "AGG")
        cand = self._cand(g, GUIDE, flank)
        annotate_specificity([cand], g, PipelineConfig())
        assert cand.specificity_score == 100.0 and cand.n_offtargets == 0

    def test_perfect_duplicate_halves_score(self):
        g, flank = genome_with(GUIDE, "AGG")
        g2 = g + rnd(200, 8) + GUIDE + "TGG" + rnd(200, 9)
        cand = self._cand(g2, GUIDE, flank)
        cand2 = [c for c in find_pam_sites(g2) if c.protospacer == GUIDE and c.strand == "+"]
        assert len(cand2) == 2
        annotate_specificity([cand2[0]], g2, PipelineConfig())
        assert cand2[0].specificity_score == pytest.approx(50.0)

    def test_seed_mismatch_copy_hurts_less_than_distal_copy(self):
        base, _ = genome_with(GUIDE, "AGG")
        seed_copy = base + rnd(100, 3) + mutate(GUIDE, [19]) + "AGG" + rnd(100, 4)
        distal_copy = base + rnd(100, 3) + mutate(GUIDE, [0]) + "AGG" + rnd(100, 4)
        scores = []
        for g in (seed_copy, distal_copy):
            cand = [c for c in find_pam_sites(g) if c.protospacer == GUIDE][0]
            annotate_specificity([cand], g, PipelineConfig())
            scores.append(cand.specificity_score)
        assert scores[0] > scores[1]      # seed mismatch = safer guide

    def test_nag_copy_hurts_less_than_ngg_copy(self):
        base, _ = genome_with(GUIDE, "AGG")
        out = []
        for pam in ("AGG", "AAG"):
            g = base + rnd(100, 3) + GUIDE + pam + rnd(100, 4)
            cand = [c for c in find_pam_sites(g) if c.protospacer == GUIDE][0]
            annotate_specificity([cand], g, PipelineConfig())
            out.append(cand.specificity_score)
        assert out[1] > out[0]
        assert PAM_FACTORS["NAG"] < PAM_FACTORS["NGG"]


class TestRankingIntegration:
    def test_rank_uses_genome_wide_score(self):
        a = PAMCandidate(0, "+", "GACGTACGTACGTACGTACG", "AGG", 17, 55.0, False,
                         specificity_score=40.0)
        b = PAMCandidate(50, "+", "GCCGTACGTACGTACGTACG", "AGG", 67, 55.0, False,
                         specificity_score=95.0)
        ranked = rank_sgRNAs([a, b], "G" * 300, PipelineConfig())
        assert ranked[0] is b

    def test_design_sgrnas_annotates_and_avoids_repeated_guides(self):
        genome = rnd(20000, 5)
        cs, ce = 8000, 12000
        up, cl, dn = genome[cs - 500:cs], genome[cs:ce], genome[ce:ce + 500]
        cfg = PipelineConfig(top_n_sgRNAs=200)
        sg = design_sgRNAs(cl, up, dn, genome, cfg, upstream_start=cs - 500)
        assert all(c.specificity_score is not None for c in sg["left"] + sg["right"])
        best_left = sg["left"][0]
        # plant a perfect copy of the top-ranked left guide far away
        planted = genome[:15000] + best_left.protospacer + "TGG" + genome[15023:]
        sg2 = design_sgRNAs(planted[cs:ce], planted[cs - 500:cs], planted[ce:ce + 500],
                            planted, cfg, upstream_start=cs - 500)
        again = [c for c in sg2["left"] if c.protospacer == best_left.protospacer]
        assert again and again[0].specificity_score <= 50.0 + 1e-6
        assert sg2["left"][0].protospacer != best_left.protospacer

    def test_legacy_model_still_works(self):
        genome = rnd(6000, 6)
        cfg = PipelineConfig(specificity_model="legacy")
        sg = design_sgRNAs(genome[2000:4000], genome[1500:2000], genome[4000:4500],
                           genome, cfg, upstream_start=1500)
        assert sg["left"] and all(c.specificity_score is None for c in sg["left"])
