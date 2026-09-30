"""Tests for joint left/right sgRNA pair ranking."""
import dataclasses
import random

import pytest

from tar_crispr.config import PAMCandidate, PipelineConfig
from tar_crispr.pair_ranking import (
    PairCandidate, rank_pairs, best_pair, _apply_arms,
    FLANK_PENALTY_PER_BP, ARM_CUTS_INTO_BGC_PENALTY,
)


def rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def cand(genome, cut, spec=100.0, strand="+", proto=None, gc=50.0, polyt=False):
    ps = proto or genome[cut - 17:cut + 3]
    return PAMCandidate(position=cut - 17, strand=strand, protospacer=ps, pam="AGG",
                        cut_position=cut, gc_percent=gc, polyt_flag=polyt,
                        specificity_score=spec)


@pytest.fixture
def setup():
    genome = rnd(8000, 5)
    cfg = PipelineConfig(mode="in-vitro")
    return genome, cfg, 2000, 6000       # genome, config, cluster_start, cluster_end


class TestWeakestLink:
    def test_pair_uses_minimum_specificity(self, setup):
        genome, cfg, cs, ce = setup
        p = PairCandidate(cand(genome, 1900, 95), cand(genome, 6100, 60), 4200)
        assert p.specificity == 60

    def test_weak_guide_loses_to_two_good_ones(self, setup):
        genome, cfg, cs, ce = setup
        sg = {"left": [cand(genome, 1900, 99), cand(genome, 1990, 98)],
              "right": [cand(genome, 6100, 40), cand(genome, 6150, 90)]}
        best = best_pair(rank_pairs(sg, genome, cs, ce, cfg))
        assert best.right.cut_position == 6150          # 90 beats 40 despite longer flank

    def test_equal_specificity_prefers_shorter_flank(self, setup):
        genome, cfg, cs, ce = setup
        g5 = lambda cut: "G" + genome[cut - 17:cut + 2]      # all guides start with G (T7)
        sg = {"left": [cand(genome, 1600, 90, proto=g5(1600)),
                       cand(genome, 1950, 90, proto=g5(1950))],
              "right": [cand(genome, 6050, 90, proto=g5(6050))]}
        best = best_pair(rank_pairs(sg, genome, cs, ce, cfg))
        assert best.left.cut_position == 1950


class TestFilters:
    def test_pairs_not_flanking_cluster_are_dropped(self, setup):
        genome, cfg, cs, ce = setup
        sg = {"left": [cand(genome, 2100, 99), cand(genome, 1900, 80)],   # 2100 cuts inside
              "right": [cand(genome, 6100, 80)]}
        pairs = rank_pairs(sg, genome, cs, ce, cfg)
        assert [p.left.cut_position for p in pairs] == [1900]

    def test_empty_side_returns_empty(self, setup):
        genome, cfg, cs, ce = setup
        assert rank_pairs({"left": [], "right": [cand(genome, 6100)]}, genome, cs, ce, cfg) == []

    def test_internal_cut_excludes_pair_but_keeps_it_listed_last(self, setup):
        genome, cfg, cs, ce = setup
        bad = cand(genome, 1900, 99)
        g = genome[:4000] + bad.protospacer + "TGG" + genome[4023:]
        sg = {"left": [bad, cand(g, 1800, 70)], "right": [cand(g, 6100, 80)]}
        pairs = rank_pairs(sg, g, cs, ce, cfg)
        assert pairs[0].left.cut_position == 1800 and not pairs[0].excluded
        assert pairs[-1].excluded and any("EXCLUDED" in n for n in pairs[-1].notes)
        assert best_pair(pairs) is pairs[0]

    def test_best_pair_none_when_all_excluded(self, setup):
        genome, cfg, cs, ce = setup
        bad = cand(genome, 1900, 99)
        g = genome[:4000] + bad.protospacer + "TGG" + genome[4023:]
        pairs = rank_pairs({"left": [bad], "right": [cand(g, 6100, 80)]}, g, cs, ce, cfg)
        assert best_pair(pairs) is None and pairs[0].excluded


class TestComponents:
    def test_flank_penalty_value(self, setup):
        genome, cfg, cs, ce = setup
        p = rank_pairs({"left": [cand(genome, 1900)], "right": [cand(genome, 6100)]},
                       genome, cs, ce, cfg)[0]
        assert p.components["flank"] == pytest.approx(-200 * FLANK_PENALTY_PER_BP)

    def test_in_vivo_polyt_penalised_in_vitro_not(self, setup):
        genome, cfg, cs, ce = setup
        sg = {"left": [cand(genome, 1900, polyt=True, proto="G" + genome[1884:1903])],
              "right": [cand(genome, 6100, proto="G" + genome[6084:6103])]}
        vitro = rank_pairs(sg, genome, cs, ce, PipelineConfig(mode="in-vitro"))[0]
        vivo = rank_pairs(sg, genome, cs, ce, PipelineConfig(mode="in-vivo"),
                          vector_seq=None)[0]
        assert vitro.components["guide_quality"] == 0
        assert vivo.components["guide_quality"] < 0

    def test_ranks_are_consecutive(self, setup):
        genome, cfg, cs, ce = setup
        sg = {"left": [cand(genome, 1900), cand(genome, 1950)],
              "right": [cand(genome, 6100), cand(genome, 6150)]}
        assert [p.rank for p in rank_pairs(sg, genome, cs, ce, cfg)] == [1, 2, 3, 4]


class TestArmsStage:
    def test_arms_evaluated_only_for_top_n(self, setup):
        genome, cfg, cs, ce = setup
        cfg.n_pairs = 2
        vec = rnd(3000, 77)
        sg = {"left": [cand(genome, 1900 + 10 * i, 90 - i) for i in range(3)],
              "right": [cand(genome, 6100 + 10 * i, 90 - i) for i in range(2)]}
        pairs = rank_pairs(sg, genome, cs, ce, cfg, vector_seq=vec)
        assert sum(p.arms_evaluated for p in pairs) == 2
        assert all("arms" in p.components for p in pairs if p.arms_evaluated)
        assert all(any("not evaluated" in n for n in p.notes)
                   for p in pairs if not p.arms_evaluated and not p.excluded)

    def test_no_vector_means_no_arm_stage(self, setup):
        genome, cfg, cs, ce = setup
        pairs = rank_pairs({"left": [cand(genome, 1900)], "right": [cand(genome, 6100)]},
                           genome, cs, ce, cfg)
        assert not pairs[0].arms_evaluated and "arms" not in pairs[0].components

    @pytest.mark.parametrize("side", ["left", "right"])
    def test_arm_shift_into_bgc_is_penalised_on_both_sides(self, setup, monkeypatch, side):
        genome, cfg, cs, ce = setup
        from tar_crispr import pair_ranking as pr
        from tar_crispr.config import HomologyArm

        def fake_arms(fragment, g, v, c):
            n = len(fragment)
            mk = lambda s, e, end: HomologyArm(fragment[s:e], end, e - s, 50.0, s, e, True, True, [])
            left = mk(0, 50, "left")
            right = mk(n - 50, n, "right")
            # shift the chosen side by 300 bp, which reaches past the 100 bp flank
            if side == "left":
                left = mk(300, 350, "left")
            else:
                right = mk(n - 350, n - 300, "right")
            return {"left": left, "right": right}

        monkeypatch.setattr(pr, "design_homology_arms", fake_arms)
        l, r = cand(genome, 1900), cand(genome, 6100)
        pair = PairCandidate(l, r, r.cut_position - l.cut_position)
        pair.components = {"specificity": 100.0}
        pr._apply_arms(pair, genome, "ACGT" * 10, cs, cfg, ce)
        assert any(f"{side} arm shift removes" in n for n in pair.notes)
        assert pair.components["arms"] <= -ARM_CUTS_INTO_BGC_PENALTY


class TestLazyArmEvaluation:
    """The pair chosen must not depend on an arbitrary evaluation cut-off."""

    def _fake(self, monkeypatch, bad_len):
        from tar_crispr import pair_ranking as pr
        from tar_crispr.config import HomologyArm

        def fake_arms(fragment, g, v, c):
            issues = ["GC% too high", "structure", "repeat"] if len(fragment) == bad_len else []
            mk = lambda s, e, end: HomologyArm(fragment[s:e], end, e - s, 50.0, s, e,
                                               True, True, list(issues))
            n = len(fragment)
            return {"left": mk(0, 50, "left"), "right": mk(n - 50, n, "right")}
        monkeypatch.setattr(pr, "design_homology_arms", fake_arms)

    def test_better_pair_with_lower_preliminary_score_can_win(self, setup, monkeypatch):
        genome, cfg, cs, ce = setup
        cfg.n_pairs = 1                      # old behaviour: only the top pair was checked
        g5 = lambda cut: "G" + genome[cut - 17:cut + 2]
        a, b = cand(genome, 1900, 100, proto=g5(1900)), cand(genome, 1800, 96, proto=g5(1800))
        r = cand(genome, 6100, 100, proto=g5(6100))
        self._fake(monkeypatch, bad_len=6100 - 1900)        # pair A: three arm problems
        pairs = rank_pairs({"left": [a, b], "right": [r]}, genome, cs, ce, cfg,
                           vector_seq="ACGT" * 10)
        best = best_pair(pairs)
        assert best.left.cut_position == 1800               # B wins once A's arms are penalised
        assert all(p.arms_evaluated for p in pairs)

    def test_stops_once_unevaluated_pairs_cannot_win(self, setup, monkeypatch):
        genome, cfg, cs, ce = setup
        cfg.n_pairs = 1
        g5 = lambda cut: "G" + genome[cut - 17:cut + 2]
        lefts = [cand(genome, 1900, 100, proto=g5(1900)),
                 cand(genome, 1850, 60, proto=g5(1850)),
                 cand(genome, 1800, 55, proto=g5(1800))]
        r = cand(genome, 6100, 100, proto=g5(6100))
        self._fake(monkeypatch, bad_len=-1)                 # no arm problems anywhere
        pairs = rank_pairs({"left": lefts, "right": [r]}, genome, cs, ce, cfg,
                           vector_seq="ACGT" * 10)
        assert [p.arms_evaluated for p in pairs] == [True, False, False]
        assert best_pair(pairs).left.cut_position == 1900

    def test_min_n_pairs_still_evaluated_for_alternatives(self, setup, monkeypatch):
        genome, cfg, cs, ce = setup
        cfg.n_pairs = 3
        g5 = lambda cut: "G" + genome[cut - 17:cut + 2]
        lefts = [cand(genome, 1900 - 20 * i, 100 - 10 * i, proto=g5(1900 - 20 * i)) for i in range(5)]
        self._fake(monkeypatch, bad_len=-1)
        pairs = rank_pairs({"left": lefts, "right": [cand(genome, 6100, 100, proto=g5(6100))]},
                           genome, cs, ce, cfg, vector_seq="ACGT" * 10)
        assert sum(p.arms_evaluated for p in pairs) == 3
