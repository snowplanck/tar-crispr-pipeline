"""Genome-wide, position-weighted sgRNA specificity scoring.

The old ranking counted look-alike sites with a uniform mismatch tolerance and
ignored the PAM. Cas9 tolerates mismatches far more at the PAM-distal end than
in the PAM-proximal seed, and needs an NGG (weakly NAG/NGA) PAM, so this module
scores each off-target site the way the MIT server (Hsu et al. 2013) does:

    f = prod(1 - W[pos])                    over the mismatched positions
        * 1 / (((19 - d) / 19) * 4 + 1)     d = mean gap between mismatches
        * 1 / n_mismatches**2

and combines the sites of one guide as

    specificity = 100 / (1 + sum(f_i * pam_factor_i))          (0-100)

so a guide with no off-targets scores 100 and a guide with one perfect extra
copy scores 50.

IMPORTANT: ``MIT_WEIGHTS`` and ``PAM_FACTORS`` were written from memory of the
published MIT/CRISPOR implementation, not copied from a file. Check them against
Hsu et al. 2013 (Nat Biotechnol 31:827) or the CRISPOR source before citing
absolute scores in a paper. The *ordering* behaviour (seed mismatches matter
more than distal ones) is what the ranking relies on, and that is tested.

How the genome is searched
--------------------------
Only the guides being ranked matter, so instead of indexing the whole genome
(millions of sites) the genome is walked once per strand and every PAM-adjacent
site is looked up in a small dictionary of the wanted seeds (exact seed and all
single-substitution seed variants). Memory stays flat; a 7 Mb genome takes a
few seconds. Sites with two or more seed mismatches are not enumerated; their
contribution to the score is negligible.

Coordinates: 0-based; a cut is a boundary index between two bases on the
forward strand, the same convention as ``PAMCandidate.cut_position``.
"""
from __future__ import annotations

import re
from itertools import product
from typing import Iterable, Optional

from tar_crispr.config import CutSite, PAMCandidate, PipelineConfig

# Hsu et al. 2013 mismatch tolerance weights. Index 0 = PAM-distal (5' end of
# the protospacer), index 19 = PAM-proximal.
MIT_WEIGHTS = (0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079,
               0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804,
               0.685, 0.583)

# Relative cleavage of non-canonical PAMs vs NGG. Heuristic (NAG is reported at
# roughly a fifth of NGG; NGA lower). Exposed so it can be tuned.
PAM_FACTORS = {"NGG": 1.0, "NAG": 0.2, "NGA": 0.1}

_PAM_RE = re.compile(r"(?=([ACGT](?:GG|AG|GA)))")
_COMP = str.maketrans("ACGT", "TGCA")
_VALID = set("ACGT")


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _pam_class(pam: str) -> str:
    return {"GG": "NGG", "AG": "NAG", "GA": "NGA"}[pam[1:]]


def mismatch_positions(guide: str, site: str) -> tuple:
    """Positions (0 = PAM-distal) where *site* differs from *guide*."""
    return tuple(i for i, (a, b) in enumerate(zip(guide, site)) if a != b)


def offtarget_hit_score(positions: Iterable[int]) -> float:
    """MIT-style score in [0, 1] for one site with mismatches at *positions*.

    No mismatches (a perfect extra copy) scores 1.0. Mismatches near the PAM
    lower the score much more than distal ones.
    """
    pos = sorted(positions)
    n = len(pos)
    if n == 0:
        return 1.0
    s1 = 1.0
    for p in pos:
        s1 *= 1.0 - MIT_WEIGHTS[p]
    if n > 1:
        d = (pos[-1] - pos[0]) / (n - 1)
        s2 = 1.0 / (((19.0 - d) / 19.0) * 4.0 + 1.0)
    else:
        s2 = 1.0
    return s1 * s2 / (n * n)


def specificity_from_hits(weighted_scores: Iterable[float]) -> float:
    """Combine per-site scores (already multiplied by the PAM factor)."""
    return 100.0 / (1.0 + sum(weighted_scores))


def scan_offtargets(protospacers: Iterable[str],
                    genome_seq: str,
                    max_mismatches: int = 3,
                    seed_len: int = 12) -> dict:
    """Find every PAM-adjacent near-match of each protospacer in *genome_seq*.

    Returns ``{protospacer: [CutSite, ...]}``. Includes the intended site (the
    caller removes it), so a unique guide has exactly one entry.
    """
    guides = sorted({g.upper() for g in protospacers if g})
    out: dict = {g: [] for g in guides}
    if not guides:
        return out
    plen = len(guides[0])

    wanted: dict = {}
    for g in guides:
        seed = g[-seed_len:]
        wanted.setdefault(seed, []).append(g)
        if max_mismatches >= 1:
            for pos, base in product(range(seed_len), "ACGT"):
                if seed[pos] != base:
                    v = seed[:pos] + base + seed[pos + 1:]
                    wanted.setdefault(v, []).append(g)

    fwd = genome_seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")
    L = len(fwd)
    for strand, text in (("+", fwd), ("-", _rc(fwd))):
        get = wanted.get
        for m in _PAM_RE.finditer(text):
            i = m.start()
            if i < plen:
                continue
            gl = get(text[i - seed_len:i])
            if gl is None:
                continue
            site = text[i - plen:i]
            if not _VALID.issuperset(site):
                continue
            pam = text[i:i + 3]
            cut = i - 3 if strand == "+" else L - (i - 3)
            for g in gl:
                pos = mismatch_positions(g, site)
                if len(pos) > max_mismatches:
                    continue
                out[g].append(CutSite(
                    cut_position=cut, strand=strand, pam=pam,
                    pam_class=_pam_class(pam), mismatches=len(pos),
                    seed_mismatches=sum(1 for p in pos if p >= plen - seed_len),
                    severity="high" if not pos or pos[-1] < plen - seed_len else "medium",
                    locus="genome", site=site,
                    score=offtarget_hit_score(pos), mismatch_positions=pos))
    return out


def annotate_specificity(candidates: list[PAMCandidate],
                         genome_seq: str,
                         config: PipelineConfig) -> None:
    """Set ``specificity_score`` / ``n_offtargets`` on every candidate in place.

    ``cut_position`` must be in ``genome_seq`` coordinates: the intended site is
    recognised by (strand, cut) and excluded from the score.
    """
    if not candidates:
        return
    if any(len(c.protospacer) != len(MIT_WEIGHTS) for c in candidates):
        # The positional weights only exist for 20-nt protospacers.
        return
    scans = scan_offtargets((c.protospacer for c in candidates), genome_seq,
                            config.max_mismatches, config.seed_length)
    for c in candidates:
        hits = scans.get(c.protospacer.upper(), [])
        others = [h for h in hits
                  if not (h.strand == c.strand and h.cut_position == c.cut_position)]
        if len(others) == len(hits) and hits:
            # The intended site was not found at the candidate's own
            # coordinates: genome_seq and cut_position disagree (e.g. a local
            # vs genomic coordinate mix-up). The score would silently count the
            # guide as its own off-target, so say so.
            c.warnings.append(
                "intended site not found at its expected coordinates; "
                "specificity may be underestimated (coordinate mismatch?)")
        weighted = [h.score * PAM_FACTORS[h.pam_class] for h in others]
        c.specificity_score = round(specificity_from_hits(weighted), 2)
        c.n_offtargets = len(others)


__all__ = [
    "MIT_WEIGHTS", "PAM_FACTORS", "mismatch_positions", "offtarget_hit_score",
    "specificity_from_hits", "scan_offtargets", "annotate_specificity",
]
