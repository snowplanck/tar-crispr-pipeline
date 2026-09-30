"""Cas9 cut-site scanning and safety checks for chosen sgRNAs.

Why this exists
---------------
A guide that also cleaves *inside* the BGC fragment (or, in vivo, inside the
capture vector or the yeast genome) ruins the capture. ``check_specificity`` in
``pam_finder`` only counts protospacer look-alikes for ranking; this module asks
the sharper question: *where would Cas9 actually cut, and is that a problem in
the chosen delivery mode?*

Model (heuristic, deliberately conservative)
--------------------------------------------
A genomic site is a potential cut site when it carries a PAM and the guide
matches it with few mismatches, weighting the PAM-proximal seed (default the 12
nt next to the PAM) much more than the distal end:

* ``high``   : NGG PAM, perfect seed, total mismatches <= ``max_mismatches``.
* ``medium`` : NAG/NGA PAM with perfect seed and <= 2 mismatches, or an NGG PAM
  with exactly one seed mismatch and total <= ``max_mismatches``.

Sites are found through a dictionary keyed on the seed, so each guide costs a
handful of lookups instead of a scan of the whole region.

Coordinates follow the rest of the code base: 0-based, a cut is the boundary
index between two bases (``seq[:cut]`` | ``seq[cut:]``), always on the forward
strand of the sequence that was indexed (plus ``offset``).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from itertools import product
from typing import Iterable, Optional

from tar_crispr.config import CutSite, PAMCandidate, PipelineConfig

_COMP = str.maketrans("ACGT", "TGCA")
_PAM_RE = re.compile(r"(?=([ACGT](?:GG|AG|GA)))")
_VALID = set("ACGT")
_MODES = ("in-vitro", "in-vivo")


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _pam_class(pam: str) -> str:
    tail = pam[1:]
    return {"GG": "NGG", "AG": "NAG", "GA": "NGA"}[tail]


@dataclass
class SiteIndex:
    """Seed-keyed index of every PAM-adjacent protospacer in a sequence."""
    seed_len: int
    protospacer_len: int
    label: str = "fragment"
    seed_map: dict = field(default_factory=dict)
    n_sites: int = 0

    def add(self, other: "SiteIndex") -> None:
        for seed, entries in other.seed_map.items():
            self.seed_map.setdefault(seed, []).extend(entries)
        self.n_sites += other.n_sites


def build_site_index(seq: str,
                     offset: int = 0,
                     protospacer_len: int = 20,
                     seed_len: int = 12,
                     circular: bool = False,
                     label: str = "fragment",
                     index: Optional[SiteIndex] = None) -> SiteIndex:
    """Index all NGG/NAG/NGA-adjacent protospacers on both strands of *seq*.

    ``offset`` is added to every reported cut position (use the genomic start of
    *seq* so cuts come out in genome coordinates). For a plasmid pass
    ``circular=True`` so sites spanning the origin are found.
    """
    seq = seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")
    n = len(seq)
    idx = index if index is not None else SiteIndex(seed_len, protospacer_len, label)
    if n < protospacer_len + 3:
        return idx

    scan = seq + seq[:protospacer_len + 2] if circular else seq
    L = len(scan)
    seen: set = set()

    for strand, text in (("+", scan), ("-", _rc(scan))):
        for m in _PAM_RE.finditer(text):
            i = m.start()
            if i < protospacer_len:
                continue
            ps = text[i - protospacer_len:i]
            if not _VALID.issuperset(ps):
                continue
            pam = text[i:i + 3]
            cut_local = i - 3 if strand == "+" else L - (i - 3)
            if circular:
                start_fwd = (i - protospacer_len) if strand == "+" else (L - i)
                key = (strand, start_fwd % n)
                if key in seen:
                    continue
                seen.add(key)
                cut_local %= n
            cut = cut_local + offset
            idx.seed_map.setdefault(ps[-seed_len:], []).append(
                (ps, pam, _pam_class(pam), cut, strand))
            idx.n_sites += 1
    return idx


def build_index_from_fasta(path: str, protospacer_len: int = 20,
                           seed_len: int = 12, label: str = "yeast") -> SiteIndex:
    """Index every record of a FASTA file (e.g. the S. cerevisiae genome)."""
    from Bio import SeqIO
    idx = SiteIndex(seed_len, protospacer_len, label)
    for rec in SeqIO.parse(path, "fasta"):
        build_site_index(str(rec.seq), 0, protospacer_len, seed_len,
                         circular=False, label=label, index=idx)
    return idx


def _seed_variants(seed: str) -> Iterable[str]:
    """All sequences at Hamming distance exactly 1 from *seed*."""
    for pos, base in product(range(len(seed)), "ACGT"):
        if seed[pos] != base:
            yield seed[:pos] + base + seed[pos + 1:]


def find_cut_sites(protospacer: str,
                   index: SiteIndex,
                   max_mismatches: int = 3) -> list[CutSite]:
    """Return every potential cut site of *protospacer* in *index*."""
    ps = protospacer.upper()
    sl = index.seed_len
    seed, distal = ps[-sl:], ps[:-sl]
    hits: list[CutSite] = []

    def _scan(seed_key: str, seed_mm: int) -> None:
        for site, pam, pclass, cut, strand in index.seed_map.get(seed_key, ()):
            dmm = sum(1 for a, b in zip(site[:-sl], distal) if a != b)
            total = dmm + seed_mm
            severity = None
            if pclass == "NGG" and seed_mm == 0 and total <= max_mismatches:
                severity = "high"
            elif pclass != "NGG" and seed_mm == 0 and total <= min(2, max_mismatches):
                severity = "medium"
            elif pclass == "NGG" and seed_mm == 1 and total <= max_mismatches:
                severity = "medium"
            if severity:
                hits.append(CutSite(
                    cut_position=cut, strand=strand, pam=pam, pam_class=pclass,
                    mismatches=total, seed_mismatches=seed_mm,
                    severity=severity, locus=index.label, site=site))

    _scan(seed, 0)
    for variant in _seed_variants(seed):
        _scan(variant, 1)
    hits.sort(key=lambda h: (h.severity != "high", h.mismatches, h.cut_position))
    return hits


def sites_between(sites: list[CutSite], lo: int, hi: int) -> list[CutSite]:
    """Sites whose cut lies strictly between *lo* and *hi* (both exclusive)."""
    return [s for s in sites if lo < s.cut_position < hi]


# --------------------------------------------------------------------------
# Pair-level validation
# --------------------------------------------------------------------------
@dataclass
class PairReport:
    """Outcome of the cut-site safety checks for one left/right guide pair."""
    mode: str
    hard_issues: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    sites: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.hard_issues


def _describe(guide: str, s: CutSite, where: str) -> str:
    return (f"{guide} guide could also cut {where} at {s.cut_position} "
            f"({s.strand} strand, PAM {s.pam}, {s.mismatches} mismatch(es), "
            f"{s.seed_mismatches} in seed)")


def region_index(genome_seq: str, lo_cut: int, hi_cut: int,
                 config: PipelineConfig) -> SiteIndex:
    """Index the genomic window that will end up inside the fragment."""
    pad = config.protospacer_length + 3
    start = max(0, lo_cut - pad)
    end = min(len(genome_seq), hi_cut + pad)
    return build_site_index(genome_seq[start:end], offset=start,
                            protospacer_len=config.protospacer_length,
                            seed_len=config.seed_length, label="fragment")


def validate_guide_pair(left: PAMCandidate,
                        right: PAMCandidate,
                        genome_seq: str,
                        config: PipelineConfig,
                        vector_seq: Optional[str] = None,
                        yeast_index: Optional[SiteIndex] = None,
                        frag_index: Optional[SiteIndex] = None) -> PairReport:
    """Check both guides for cleavage inside the fragment (and, in vivo, elsewhere).

    * Fragment (always): a high-confidence site strictly between the two
      intended cuts fragments the BGC -> hard issue; medium -> warning.
    * Vector (in-vivo only): the vector is exposed to Cas9 in yeast.
    * Yeast genome (in-vivo only): requires ``yeast_index``; otherwise a
      warning says it was not evaluated.
    In ``in-vitro`` mode the vector and the host genome never meet the RNP, so
    those checks are skipped and noted.
    """
    if config.mode not in _MODES:
        raise ValueError(f"mode must be one of {_MODES}, got {config.mode!r}")
    lo, hi = left.cut_position, right.cut_position
    report = PairReport(mode=config.mode)
    if frag_index is None:
        frag_index = region_index(genome_seq, lo, hi, config)

    for name, g in (("left", left), ("right", right)):
        for s in sites_between(find_cut_sites(g.protospacer, frag_index,
                                              config.max_mismatches), lo, hi):
            report.sites.append(s)
            msg = _describe(name, s, "inside the fragment")
            (report.hard_issues if s.severity == "high" else report.warnings).append(msg)

    if config.mode == "in-vitro":
        report.notes.append(
            "in-vitro mode: Cas9 digests genomic DNA before transformation, so "
            "the capture vector and the yeast genome are not exposed; those "
            "checks were skipped.")
        return report

    # ---- in-vivo -------------------------------------------------------
    if vector_seq:
        vidx = build_site_index(vector_seq, 0, config.protospacer_length,
                                config.seed_length, circular=True, label="vector")
        for name, g in (("left", left), ("right", right)):
            for s in find_cut_sites(g.protospacer, vidx, config.max_mismatches):
                report.sites.append(s)
                msg = _describe(name, s, "in the capture vector")
                (report.hard_issues if s.severity == "high" else report.warnings).append(msg)
    else:
        report.warnings.append("in-vivo mode: no vector sequence given; vector cleavage not evaluated.")

    if yeast_index is not None:
        for name, g in (("left", left), ("right", right)):
            for s in find_cut_sites(g.protospacer, yeast_index, config.max_mismatches):
                report.sites.append(s)
                msg = _describe(name, s, "in the yeast genome")
                (report.hard_issues if s.severity == "high" else report.warnings).append(msg)
    else:
        report.warnings.append(
            "in-vivo mode: yeast genome not provided (--yeast-genome); "
            "host off-target cleavage was NOT evaluated.")
    return report


def select_valid_pair(sgRNAs: dict,
                      genome_seq: str,
                      config: PipelineConfig,
                      vector_seq: Optional[str] = None,
                      yeast_index: Optional[SiteIndex] = None):
    """Pick the best-ranked left/right pair with no hard cut-site issues.

    Pairs are tried in order of combined rank (i + j). Returns
    ``(pair_dict | None, PairReport | None, n_tried)``. When nothing passes,
    the returned report is that of the best-ranked pair (for diagnostics).
    """
    lefts, rights = sgRNAs.get("left", []), sgRNAs.get("right", [])
    if not lefts or not rights:
        return None, None, 0
    lo = min(g.cut_position for g in lefts)
    hi = max(g.cut_position for g in rights)
    frag_index = region_index(genome_seq, lo, hi, config)

    order = sorted(((i + j, i, j) for i in range(len(lefts)) for j in range(len(rights))))
    first_report = None
    tried = 0
    for _, i, j in order:
        if lefts[i].cut_position >= rights[j].cut_position:
            continue
        rep = validate_guide_pair(lefts[i], rights[j], genome_seq, config,
                                  vector_seq, yeast_index, frag_index)
        tried += 1
        first_report = first_report or rep
        if rep.ok:
            return {"left": lefts[i], "right": rights[j]}, rep, tried
    return None, first_report, tried


__all__ = [
    "SiteIndex", "PairReport", "build_site_index", "build_index_from_fasta",
    "find_cut_sites", "sites_between", "region_index",
    "validate_guide_pair", "select_valid_pair",
]
