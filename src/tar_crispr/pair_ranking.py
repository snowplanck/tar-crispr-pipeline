"""Joint ranking of left/right sgRNA pairs.

Choosing the best guide at each end independently ignores what actually gets
built: the two cuts define one fragment, and its ends define the homology arms.
This module scores complete pairs on a 0-100 scale:

    score = min(specificity_left, specificity_right)     weakest guide wins
            - guide quality penalties                    mode-aware
            - flank penalty                              extra DNA around the BGC
            - cut-site warnings                          see cut_specificity
            - homology-arm penalties                     top pairs only

A pair with a high-confidence cut inside the fragment (or, in vivo, in the
vector / yeast genome) is *excluded*: it stays in the list, at the bottom, with
the reason, so the report can explain why it was not chosen.

Homology arms are expensive (uniqueness search, secondary structure), so they
are designed lazily, best preliminary score first, until no unevaluated pair
can beat the best evaluated one (arm penalties only subtract, so the preliminary
score is an upper bound) and at least ``config.n_pairs`` pairs were evaluated.
The remaining pairs are listed after them, marked "arms not evaluated". Every
weight below is a transparent heuristic constant, not a fitted model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from tar_crispr.config import PAMCandidate, PipelineConfig
from tar_crispr.cut_specificity import (
    PairReport, SiteIndex, region_index, validate_guide_pair,
)
from tar_crispr.homology_arms import design_homology_arms

# ---- weights (points on the 0-100 scale) ---------------------------------
POLYT_PENALTY_IN_VIVO = 10.0      # Pol III terminator inside the guide
NO_5G_PENALTY_IN_VITRO = 2.0      # T7 transcription prefers a 5' G
GC_PENALTY = 3.0                  # guide GC% outside 25-75
FLANK_PENALTY_PER_BP = 0.005      # per bp of DNA beyond the BGC, capped
FLANK_PENALTY_CAP = 10.0
CUT_WARNING_PENALTY = 5.0         # per medium-confidence cut site
CUT_WARNING_CAP = 15.0
ARM_ISSUE_PENALTY = 5.0           # per arm validation issue
ARM_NONUNIQUE_PENALTY = 10.0      # arm repeats elsewhere / in the vector
ARM_SHIFT_PENALTY_PER_BP = 0.1    # bases of the fragment end given up
ARM_SHIFT_CAP = 10.0
ARM_CUTS_INTO_BGC_PENALTY = 20.0  # shifted arm gives up bases of the BGC itself


@dataclass
class PairCandidate:
    left: PAMCandidate
    right: PAMCandidate
    fragment_length: int
    score: float = 0.0
    excluded: bool = False
    components: dict = field(default_factory=dict)   # name -> points (negative = penalty)
    notes: list = field(default_factory=list)
    cut_report: Optional[PairReport] = None
    arms: Optional[dict] = None
    arms_evaluated: bool = False
    rank: int = 0

    @property
    def specificity(self) -> Optional[float]:
        vals = [g.specificity_score for g in (self.left, self.right)
                if g.specificity_score is not None]
        return min(vals) if vals else None


def _guide_penalties(g: PAMCandidate, config: PipelineConfig) -> tuple[float, list]:
    pts, notes = 0.0, []
    if config.mode == "in-vivo" and g.polyt_flag:
        pts += POLYT_PENALTY_IN_VIVO
        notes.append(f"{g.protospacer}: poly-T run (Pol III terminator)")
    if config.mode == "in-vitro" and not g.protospacer.upper().startswith("G"):
        pts += NO_5G_PENALTY_IN_VITRO
        notes.append(f"{g.protospacer}: no 5' G (prepend G for T7 transcription)")
    if g.gc_percent > 75 or g.gc_percent < 25:
        pts += GC_PENALTY
        notes.append(f"{g.protospacer}: GC {g.gc_percent}% outside 25-75%")
    return pts, notes


def _preliminary(pair: PairCandidate, cluster_start: Optional[int],
                 cluster_end: Optional[int], config: PipelineConfig,
                 genome_seq: str, vector_seq: Optional[str],
                 yeast_index: Optional[SiteIndex],
                 frag_index: SiteIndex) -> None:
    comp: dict = {}
    spec = pair.specificity
    if spec is None:
        pair.notes.append("specificity not scored (legacy model): counted as 100")
        spec = 100.0
    comp["specificity"] = spec

    pts, notes = 0.0, []
    for g in (pair.left, pair.right):
        p, n = _guide_penalties(g, config)
        pts += p
        notes += n
    comp["guide_quality"] = -pts
    pair.notes += notes

    if cluster_start is not None and cluster_end is not None:
        flank = (cluster_start - pair.left.cut_position) + (pair.right.cut_position - cluster_end)
        comp["flank"] = -min(FLANK_PENALTY_CAP, max(0, flank) * FLANK_PENALTY_PER_BP)
    else:
        comp["flank"] = 0.0

    rep = validate_guide_pair(pair.left, pair.right, genome_seq, config,
                              vector_seq, yeast_index, frag_index)
    pair.cut_report = rep
    n_warn = len(rep.warnings)
    comp["cut_warnings"] = -min(CUT_WARNING_CAP, n_warn * CUT_WARNING_PENALTY)
    if not rep.ok:
        pair.excluded = True
        pair.notes += [f"EXCLUDED: {m}" for m in rep.hard_issues]
    pair.notes += rep.warnings

    pair.components = comp
    pair.score = sum(comp.values())


def _apply_arms(pair: PairCandidate, genome_seq: str, vector_seq: str,
                cluster_start: Optional[int], config: PipelineConfig,
                cluster_end: Optional[int] = None) -> None:
    fragment = genome_seq[pair.left.cut_position:pair.right.cut_position]
    arms = design_homology_arms(fragment, genome_seq, vector_seq, config)
    pair.arms, pair.arms_evaluated = arms, True

    pen = 0.0
    shifted = {}
    for side in ("left", "right"):
        arm = arms.get(side)
        if arm is None:
            pen += ARM_NONUNIQUE_PENALTY
            pair.notes.append(f"{side} arm: no arm could be designed")
            continue
        pen += ARM_ISSUE_PENALTY * len(arm.issues)
        if not arm.uniqueness:
            pen += ARM_NONUNIQUE_PENALTY
        for issue in arm.issues:
            pair.notes.append(f"{side} arm: {issue}")
        shifted[side] = arm.start_pos if side == "left" else len(fragment) - arm.end_pos

    shift_total = sum(shifted.values())
    pen += min(ARM_SHIFT_CAP, shift_total * ARM_SHIFT_PENALTY_PER_BP)
    if shift_total:
        pair.notes.append(
            f"arms shifted inward by {shifted.get('left', 0)} bp (left) / "
            f"{shifted.get('right', 0)} bp (right): those end bases are not captured")

    # A shifted arm that reaches into the BGC drops real cluster sequence.
    if cluster_start is not None and "left" in shifted:
        loss = max(0, pair.left.cut_position + shifted["left"] - cluster_start)
        if loss:
            pen += ARM_CUTS_INTO_BGC_PENALTY
            pair.notes.append(f"left arm shift removes {loss} bp of the BGC itself")
    if cluster_end is not None and "right" in shifted:
        loss = max(0, cluster_end - (pair.right.cut_position - shifted["right"]))
        if loss:
            pen += ARM_CUTS_INTO_BGC_PENALTY
            pair.notes.append(f"right arm shift removes {loss} bp of the BGC itself")

    pair.components["arms"] = -pen
    pair.score = sum(pair.components.values())


ARM_EVAL_HARD_CAP = 40      # never design arms for more pairs than this


def _evaluate_arms_lazily(live: list, genome_seq: str, vector_seq: str,
                          cluster_start: Optional[int], cluster_end: Optional[int],
                          config: PipelineConfig) -> list:
    """Evaluate homology arms in order of preliminary score, stopping when safe.

    Arm penalties only subtract, so a pair's preliminary score is an upper
    bound on its final score. Pairs are therefore evaluated best-first until
    (a) at least ``config.n_pairs`` were evaluated, so the report has
    alternatives, and (b) the best evaluated score is >= the preliminary score
    of every remaining pair, which proves no unevaluated pair can win. Stopping
    earlier could pick a pair that only looks best because the better ones
    were never checked. ``ARM_EVAL_HARD_CAP`` bounds the runtime; if it is
    reached before the proof completes, a note says so.
    """
    pending = list(live)                       # already sorted by preliminary score
    evaluated: list = []
    while pending and len(evaluated) < ARM_EVAL_HARD_CAP:
        if len(evaluated) >= max(1, config.n_pairs):
            best_final = max(p.score for p in evaluated)
            if best_final >= pending[0].score:
                break
        nxt = pending.pop(0)
        _apply_arms(nxt, genome_seq, vector_seq, cluster_start, config, cluster_end)
        evaluated.append(nxt)

    evaluated.sort(key=lambda p: -p.score)
    if pending and evaluated and evaluated[0].score < pending[0].score:
        for p in evaluated[:1]:
            p.notes.append(
                f"arm evaluation stopped at {ARM_EVAL_HARD_CAP} pairs; an unevaluated pair "
                "could still score higher (raise --top-n less, or narrow the flank window)")
    for p in pending:
        p.notes.append("homology arms not evaluated for this pair "
                       "(its best possible score is below the selected pair)")
    return evaluated + pending


def rank_pairs(sgRNAs: dict,
               genome_seq: str,
               cluster_start: Optional[int],
               cluster_end: Optional[int],
               config: PipelineConfig,
               vector_seq: Optional[str] = None,
               yeast_index: Optional[SiteIndex] = None) -> list[PairCandidate]:
    """Return all valid left/right pairs, best first.

    ``sgRNAs`` is the ``{"left": [...], "right": [...]}`` dict from
    ``design_sgRNAs``. When ``cluster_start``/``cluster_end`` are given, pairs
    whose cuts do not flank the cluster are dropped. Excluded pairs (cut-site
    hard issues) are kept at the end with ``excluded=True``.
    """
    lefts, rights = sgRNAs.get("left", []), sgRNAs.get("right", [])
    if not lefts or not rights:
        return []

    lo = min(g.cut_position for g in lefts)
    hi = max(g.cut_position for g in rights)
    frag_index = region_index(genome_seq, lo, hi, config)

    pairs: list[PairCandidate] = []
    for l in lefts:
        for r in rights:
            if l.cut_position >= r.cut_position:
                continue
            if cluster_start is not None and l.cut_position > cluster_start:
                continue
            if cluster_end is not None and r.cut_position < cluster_end:
                continue
            pair = PairCandidate(left=l, right=r,
                                 fragment_length=r.cut_position - l.cut_position)
            _preliminary(pair, cluster_start, cluster_end, config, genome_seq,
                         vector_seq, yeast_index, frag_index)
            pairs.append(pair)

    live = sorted((p for p in pairs if not p.excluded), key=lambda p: -p.score)
    dead = sorted((p for p in pairs if p.excluded), key=lambda p: -p.score)

    if vector_seq:
        live = _evaluate_arms_lazily(live, genome_seq, vector_seq,
                                     cluster_start, cluster_end, config)

    ordered = live + dead
    for i, p in enumerate(ordered, 1):
        p.rank = i
    return ordered


def best_pair(pairs: list[PairCandidate]) -> Optional[PairCandidate]:
    """First non-excluded pair, or None."""
    return next((p for p in pairs if not p.excluded), None)


__all__ = ["PairCandidate", "rank_pairs", "best_pair"]
