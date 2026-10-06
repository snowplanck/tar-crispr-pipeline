"""Colony-PCR screening primers for the assembled TAR construct.

After TAR cloning in yeast most colonies carry the empty re-circularised vector
or a partial/rearranged insert, so clones are screened by PCR before anything
is sequenced. This module designs that panel on the *final construct*:

* **junction-left / junction-right** - one primer in the vector backbone and one
  inside the BGC fragment, so a band appears only when the vector is joined to
  the right end of the BGC. The empty vector gives no band.
* **marker-gene** - an amplicon inside a core biosynthetic gene (type II PKS
  KS/CLF by default), found in the GenBank annotation of the BGC.
* **integrity** - amplicons spread along the BGC, to catch deletions and
  partial captures that still have both junctions. They are placed so that the
  largest stretch of BGC without an interior amplicon is as short as possible,
  taking the marker gene into account.

Product sizes differ by at least 150 bp between any two amplicons, so the panel
can be read from a single gel or run as one multiplex. Candidate pairs come from primer3; each is
then checked by an in-silico PCR on the construct (the intended product must be
the only one for that pair), against an optional host genome (yeast), and for
cross-dimers with the primers already chosen.

Coordinates are 0-based, half-open, on the forward strand of the construct.
Tm uses primer3's ``calc_tm`` defaults (50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM
primer), the same conditions the rest of the pipeline reports.
"""
from __future__ import annotations

import csv
import re
from typing import Optional

import primer3

from tar_crispr.config import (
    PipelineConfig, ScreeningAmplicon, ScreeningDesign, ScreeningPrimer,
)

_COMP = str.maketrans("ACGTN", "TGCAN")

# Product-size slots. Each slot is 100 bp wide and consecutive slots are 150 bp apart
# (top of one slot to the bottom of the next), so any two products of one panel differ by
# at least MIN_PRODUCT_SEPARATION bp and resolve as separate bands on a gel. (Slots 50 bp
# apart, as an earlier version had, let two products land 92 bp apart.)
MIN_PRODUCT_SEPARATION = 150
SLOT_WIDTH = 100


def size_slot(i: int) -> tuple:
    """The i-th product-size slot: (300, 400), (550, 650), (800, 900), ..."""
    lo = 300 + (SLOT_WIDTH + MIN_PRODUCT_SEPARATION) * i
    return (lo, lo + SLOT_WIDTH)


SIZE_LADDER = [size_slot(i) for i in range(8)]

CROSS_DIMER_WARN_DG = -9000.0     # cal/mol (primer3 heterodimer dG at 37 C)
MAX_MISMATCH_SITE = 3             # a site with <= this many mismatches can prime
ANCHOR = 10                       # exact 3'-end match required to prime
NUM_CANDIDATES = 20

# (name, overrides applied to the strict primer3 settings)
_LEVELS = [
    ("strict", {}),
    ("relaxed", {"PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0, "PRIMER_MIN_GC": 40.0,
                 "PRIMER_MAX_GC": 80.0, "PRIMER_MAX_END_GC": 4, "PRIMER_MAX_SIZE": 28}),
    ("very relaxed", {"PRIMER_MIN_TM": 55.0, "PRIMER_MAX_TM": 65.0, "PRIMER_MIN_GC": 35.0,
                      "PRIMER_MAX_GC": 85.0, "PRIMER_MAX_END_GC": 5, "PRIMER_MAX_SIZE": 30,
                      "PRIMER_GC_CLAMP": 0}),
]


def _clean(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper()


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _gc(seq: str) -> float:
    return round(100.0 * sum(1 for b in seq if b in "GC") / len(seq), 1) if seq else 0.0


# --------------------------------------------------------------------------
# In-silico PCR
# --------------------------------------------------------------------------
class Background:
    """A template sequence prepared once for many primer look-ups.

    Cleaning and reverse-complementing a 12 Mb genome for every candidate primer
    would dominate the run time, so the host genome is wrapped in one of these.
    """

    def __init__(self, seq: str, circular: bool):
        self.seq = _clean(seq)
        self.circular = circular
        self.length = len(self.seq)
        self.rc = _rc(self.seq)

    def sites(self, primer: str, max_mismatches: int = MAX_MISMATCH_SITE) -> list:
        return _find_sites(_clean(primer), self, max_mismatches)


def _sites_one_strand(primer: str, text: str, max_mm: int) -> list:
    """(start, end, mismatches) of primer-like matches in *text*, 3' end anchored."""
    n = len(primer)
    anchor = primer[-ANCHOR:]
    out = []
    i = text.find(anchor)
    while i != -1:
        end = i + ANCHOR
        start = end - n
        if start >= 0:
            mm = sum(1 for a, b in zip(primer, text[start:end]) if a != b)
            if mm <= max_mm:
                out.append((start, end, mm))
        i = text.find(anchor, i + 1)
    return out


def _find_sites(primer: str, bg: Background, max_mm: int) -> list:
    n, L = len(primer), bg.length
    if L < n:
        return []
    ext = n + ANCHOR if bg.circular else 0
    sites = []
    top = bg.seq + bg.seq[:ext] if bg.circular else bg.seq
    for s, e, mm in _sites_one_strand(primer, top, max_mm):
        if s < L:
            sites.append(("+", s, e, mm))
    rtext = bg.rc + bg.rc[:ext] if bg.circular else bg.rc
    for s, e, mm in _sites_one_strand(primer, rtext, max_mm):
        if s < L:
            t_start = (L - e) % L if bg.circular else L - e
            if t_start < 0:
                continue
            sites.append(("-", t_start, t_start + n, mm))
    return sites


def find_primer_sites(primer: str, seq: str, circular: bool = True,
                      max_mismatches: int = MAX_MISMATCH_SITE) -> list:
    """Every place *primer* can prime on either strand of *seq*.

    Returns ``(strand, start, end, mismatches)`` on the forward strand of *seq*;
    ``"+"`` means the primer is identical to the top strand (a forward primer),
    ``"-"`` means its reverse complement is (a reverse primer). A site must
    match exactly at the 3' end (last 10 nt) and have at most
    ``max_mismatches`` mismatches overall. For a circular *seq* a site may span
    the origin (then ``end`` exceeds ``len(seq)``).
    """
    return _find_sites(_clean(primer), Background(seq, circular), max_mismatches)


def _pair_products(fname: str, rname: str, f_sites: list, r_sites: list,
                   L: int, circular: bool, max_size: int) -> list:
    out = []
    for fs_strand, fs, fe, fmm in f_sites:
        if fs_strand != "+":
            continue
        for rs_strand, rs, re_, rmm in r_sites:
            if rs_strand != "-":
                continue
            if re_ > fs:
                size = re_ - fs
            elif circular:
                size = re_ + L - fs
            else:
                continue
            if size < max(fe - fs, re_ - rs) or size > max_size:
                continue
            out.append({"forward": fname, "reverse": rname, "start": fs, "end": fs + size,
                        "size": size, "mismatches": fmm + rmm})
    return out


def predict_amplicons(primers: dict, seq: str, circular: bool = True,
                      max_size: int = 4000,
                      max_mismatches: int = MAX_MISMATCH_SITE) -> list:
    """In-silico PCR: every product any forward/reverse primer combination gives.

    ``primers`` maps name -> sequence. Returns dicts with ``forward``,
    ``reverse``, ``start``, ``end`` (half-open; ``end`` may exceed ``len(seq)``
    when the product spans the origin of a circular *seq*), ``size`` and the
    total ``mismatches`` of the two sites.
    """
    bg = Background(seq, circular)
    sites = {name: bg.sites(p, max_mismatches) for name, p in primers.items()}
    out = []
    for fn, fsites in sites.items():
        for rn, rsites in sites.items():
            out += _pair_products(fn, rn, fsites, rsites, bg.length, circular, max_size)
    out.sort(key=lambda d: (d["size"], d["forward"], d["reverse"]))
    return out


# --------------------------------------------------------------------------
# Marker gene (from the BGC annotation)
# --------------------------------------------------------------------------
def _feature_text(feature) -> str:
    """Searchable annotation text of a feature (original case; lower() it to search)."""
    parts = []
    for key in ("product", "gene", "gene_functions", "sec_met_domain", "note", "function"):
        parts.extend(str(v) for v in feature.qualifiers.get(key, []))
    return re.sub(r"\s+", " ", " ".join(parts))


def _kw_regex(keyword: str):
    """Regex matching *keyword* ignoring case and treating ``-``, ``_`` and spaces alike.

    antiSMASH writes the same domain as ``Chain_length_factor``, ``chain-length-factor``
    or ``chain length factor`` depending on the field, so all three must match. Matching
    on the original text (rather than a normalised copy) lets the report quote the
    annotation exactly as it appears.
    """
    tokens = [t for t in re.split(r"[\s_\-]+", keyword.strip()) if t]
    if not tokens:
        return None
    return re.compile(r"[\s_\-]+".join(re.escape(t) for t in tokens), re.IGNORECASE)


def _cds_name(feature) -> str:
    q = feature.qualifiers
    return str((q.get("locus_tag") or q.get("gene") or q.get("protein_id") or ["CDS"])[0])


def _locate_cds(feature, bgc_record, fragment: str):
    """Place a CDS in the fragment by its nucleotide sequence (either strand).

    Returns ``((start, length, strand), None)`` or ``(None, reason)``. Searching by
    sequence needs no coordinate translation between the annotation record, the
    genome and the fragment.
    """
    try:
        gene = _clean(str(feature.extract(bgc_record.seq)))
    except Exception:
        return None, "its sequence could not be extracted from the annotation"
    if len(gene) < 90:
        return None, "it is shorter than 90 bp"
    for strand, query in (("+", gene), ("-", _rc(gene))):
        n = fragment.count(query)
        if n == 1:
            return (fragment.find(query), len(gene), strand), None
        if n > 1:
            return None, "its sequence occurs more than once in the fragment"
    return None, "its sequence is not in the fragment (annotation and fragment differ?)"


def _marker_by_tag(bgc_record, fragment: str, tag: str) -> dict:
    """The CDS named *tag* (locus_tag, gene or protein_id; case-insensitive)."""
    wanted = tag.strip().lower()
    cds = [f for f in bgc_record.features if f.type == "CDS"]
    matches = [f for f in cds
               if any(str(v).strip().lower() == wanted
                      for k in ("locus_tag", "gene", "protein_id") for v in f.qualifiers.get(k, []))]
    if not matches:
        names = [_cds_name(f) for f in sorted(cds, key=lambda f: int(f.location.start))]
        shown = ", ".join(names[:20]) + (" ..." if len(names) > 20 else "")
        raise ValueError(f"marker gene {tag!r} not found: no CDS has it as locus_tag, gene or "
                         f"protein_id. CDSs in the annotation: {shown or 'none'}")
    if len(matches) > 1:
        raise ValueError(f"marker gene {tag!r} is ambiguous ({len(matches)} CDSs); use the locus_tag")
    f = matches[0]
    loc, reason = _locate_cds(f, bgc_record, fragment)
    if loc is None:
        raise ValueError(f"marker gene {tag!r} is in the annotation but cannot be placed in the "
                         f"fragment: {reason}")
    pos, length, strand = loc
    return {"name": _cds_name(f), "product": str((f.qualifiers.get("product") or [""])[0]),
            "keyword": "--marker-gene", "evidence": f"selected by name ({tag})",
            "start": pos, "end": pos + length, "strand": strand, "n_matches": 1,
            "alternatives": []}


def find_marker_gene(bgc_record, fragment_seq: str, keywords: list,
                     locus_tag: Optional[str] = None) -> Optional[dict]:
    """Locate the marker gene of the BGC annotation inside the fragment.

    With ``locus_tag`` the named CDS is used (a ``ValueError`` explains why if it is not
    in the annotation or cannot be placed in the fragment). Otherwise keywords are tried
    in priority order against each CDS's product/gene/gene_functions/sec_met_domain/note
    text, ignoring case and treating ``-``, ``_`` and spaces alike; among CDSs matching
    the first productive keyword, the one nearest the fragment centre is used and the
    others are returned under ``alternatives``. Returns ``None`` if nothing matches.
    """
    if locus_tag:
        if bgc_record is None or not getattr(bgc_record, "features", None):
            raise ValueError("--marker-gene needs an annotated BGC (give --genbank)")
        return _marker_by_tag(bgc_record, _clean(fragment_seq), locus_tag)
    if bgc_record is None or not getattr(bgc_record, "features", None):
        return None
    fragment = _clean(fragment_seq)
    centre = len(fragment) / 2
    cds = [f for f in bgc_record.features if f.type == "CDS"]
    for kw in keywords:
        rx = _kw_regex(kw)
        kw = kw.lower()
        hits = []
        for f in cds:
            text = _feature_text(f)
            m_kw = rx.search(text) if rx else None
            if m_kw is None:
                continue
            loc, _ = _locate_cds(f, bgc_record, fragment)
            if loc is None:
                continue
            pos, length, strand = loc
            at, kw_len = m_kw.start(), m_kw.end() - m_kw.start()
            hits.append({"name": _cds_name(f),
                         "product": str((f.qualifiers.get("product") or [""])[0]),
                         "keyword": kw, "evidence": text[max(0, at - 30): at + kw_len + 50].strip(),
                         "start": pos, "end": pos + length, "strand": strand})
        if hits:
            hits.sort(key=lambda h: abs((h["start"] + h["end"]) / 2 - centre))
            best = hits[0]
            best["n_matches"] = len(hits)
            best["alternatives"] = [{k: h[k] for k in ("name", "product", "start", "end")}
                                    for h in hits[1:]]
            return best
    return None


# --------------------------------------------------------------------------
# Primer3 wrapper
# --------------------------------------------------------------------------
def _p3_globals(config: PipelineConfig, size_range: tuple, overrides: dict) -> dict:
    g = {
        "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 26,
        "PRIMER_OPT_TM": (config.primer_tm_min + config.primer_tm_max) / 2,
        "PRIMER_MIN_TM": config.primer_tm_min, "PRIMER_MAX_TM": config.primer_tm_max,
        "PRIMER_MIN_GC": 45.0, "PRIMER_MAX_GC": 75.0,
        "PRIMER_MAX_END_GC": 3, "PRIMER_GC_CLAMP": 1, "PRIMER_MAX_POLY_X": 4,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(size_range)],
        "PRIMER_NUM_RETURN": NUM_CANDIDATES,
        # Same conditions as primer3.calc_tm's defaults (50 mM monovalent, 1.5 mM Mg2+,
        # 0.6 mM dNTP, 50 nM primer), which the rest of the pipeline uses to report Tm.
        "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_SALT_DIVALENT": 1.5,
        "PRIMER_DNTP_CONC": 0.6, "PRIMER_DNA_CONC": 50.0,
    }
    g.update(overrides)
    return g


def _candidates(construct: str, win_start: int, win_end: int,
                left_region: tuple, right_region: tuple,
                config: PipelineConfig, size_range: tuple, overrides: dict) -> list:
    """primer3 pairs inside construct[win_start:win_end]; regions in construct coords."""
    win_start, win_end = max(0, win_start), min(len(construct), win_end)
    tmpl = construct[win_start:win_end]
    ls, le = max(left_region[0], win_start), min(left_region[1], win_end)
    rs, re_ = max(right_region[0], win_start), min(right_region[1], win_end)
    if le - ls < 18 or re_ - rs < 18 or len(tmpl) < size_range[0]:
        return []
    seq_args = {
        "SEQUENCE_ID": "screen", "SEQUENCE_TEMPLATE": tmpl,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [[ls - win_start, le - ls,
                                                 rs - win_start, re_ - rs]],
    }
    try:
        res = primer3.bindings.design_primers(seq_args, _p3_globals(config, size_range, overrides))
    except Exception:
        return []
    out = []
    for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        lpos, rpos = res[f"PRIMER_LEFT_{i}"], res[f"PRIMER_RIGHT_{i}"]
        fseq = res[f"PRIMER_LEFT_{i}_SEQUENCE"].upper()
        rseq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"].upper()
        f_start = win_start + lpos[0]
        r_end = win_start + rpos[0] + 1              # primer3 gives the 5' end of the reverse primer
        r_start = r_end - rpos[1]
        out.append({"f": fseq, "r": rseq, "f_start": f_start, "f_end": f_start + lpos[1],
                    "r_start": r_start, "r_end": r_end,
                    "size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                    "penalty": res[f"PRIMER_PAIR_{i}_PENALTY"]})
    return out


def _tm(seq: str) -> float:
    return round(primer3.calc_tm(seq), 1)


def _dimer_dg(a: str, b: str) -> float:
    try:
        return primer3.calc_heterodimer(a, b, temp_c=37).dg
    except Exception:
        return 0.0


def _evaluate(c: dict, cons: Background, chosen: list, host: Optional[Background],
              fwd_name: str, rev_name: str) -> tuple:
    """Return (issues, warnings, score) for one primer3 candidate pair."""
    issues, warnings = [], []
    sites = {fwd_name: cons.sites(c["f"]), rev_name: cons.sites(c["r"])}

    prods = _pair_products(fwd_name, rev_name, sites[fwd_name], sites[rev_name],
                           cons.length, True, 4000)
    intended = [p for p in prods if p["size"] == c["size"] and p["start"] == c["f_start"]]
    unexpected = [p for p in prods if p not in intended]
    if not intended:
        issues.append("intended product not reproduced by in-silico PCR")
    if unexpected:
        issues.append("unintended product(s) on the construct: "
                      + ", ".join(f"{p['size']} bp ({p['mismatches']} mm)" for p in unexpected))

    intended_site = {fwd_name: ("+", c["f_start"]), rev_name: ("-", c["r_start"])}
    for name, ss in sites.items():
        extra = [x for x in ss if not (x[3] == 0 and (x[0], x[1]) == intended_site[name])]
        if extra:
            warnings.append(f"{name}: {len(extra)} additional binding site(s) on the construct "
                            f"(<= {MAX_MISMATCH_SITE} mismatches, exact 3' 10-mer)")

    if host is not None:
        hs = {fwd_name: host.sites(c["f"]), rev_name: host.sites(c["r"])}
        for name, ss in hs.items():
            if ss:
                warnings.append(f"{name}: {len(ss)} site(s) in the host genome")
        hp = _pair_products(fwd_name, rev_name, hs[fwd_name], hs[rev_name],
                            host.length, False, 4000)
        if hp:
            warnings.append("predicted host-genome product(s): "
                            + ", ".join(f"{p['size']} bp" for p in hp[:3]))

    worst = 0.0
    for other in chosen:
        for a in (c["f"], c["r"]):
            for b in (other.forward.sequence, other.reverse.sequence):
                worst = min(worst, _dimer_dg(a, b))
    if worst < CROSS_DIMER_WARN_DG:
        warnings.append(f"cross-dimer with another screening primer (dG {worst / 1000:.1f} kcal/mol)")

    return issues, warnings, (len(issues), len(warnings), c["penalty"])


def _design_one(name: str, role: str, target: str, construct: str, cons: Background,
                win: tuple, left_region: tuple, right_region: tuple, size_range: tuple,
                config: PipelineConfig, chosen: list,
                host: Optional[Background]) -> Optional[ScreeningAmplicon]:
    fwd_name, rev_name = f"{name}_F", f"{name}_R"
    best = None
    for level, overrides in _LEVELS:
        cands = _candidates(construct, win[0], win[1], left_region, right_region,
                            config, size_range, overrides)
        scored = []
        for c in cands:                      # primer3 order = best penalty first
            issues, warnings, score = _evaluate(c, cons, chosen, host, fwd_name, rev_name)
            scored.append((score, c, issues, warnings))
            if score[0] == 0 and score[1] == 0:
                break                        # first clean pair wins; no need to check the rest
        if not scored:
            continue
        scored.sort(key=lambda t: t[0])
        cand_best = scored[0]
        if best is None or cand_best[0][0] < best[0][0][0]:
            best = (cand_best, level)
        if cand_best[0][0] == 0:         # no hard issues: stop relaxing
            best = (cand_best, level)
            break
    if best is None:
        return None
    (score, c, issues, warnings), level = best
    if level != "strict":
        warnings = warnings + [f"primers needed the '{level}' design settings "
                               "(wider Tm / GC window)"]
    fp = ScreeningPrimer(fwd_name, c["f"], "+", c["f_start"], c["f_end"],
                         _tm(c["f"]), _gc(c["f"]), len(c["f"]))
    rp = ScreeningPrimer(rev_name, c["r"], "-", c["r_start"], c["r_end"],
                         _tm(c["r"]), _gc(c["r"]), len(c["r"]))
    return ScreeningAmplicon(name=name, role=role, target=target, forward=fp, reverse=rp,
                             product_size=c["size"], start=c["f_start"], end=c["r_end"],
                             level=level, issues=issues, warnings=warnings)


# --------------------------------------------------------------------------
# Public entry point
# --------------------------------------------------------------------------
def _allocate_centres(frag_len: int, n: int, fixed: Optional[list] = None,
                      margin: int = 700) -> list:
    """Centres for *n* integrity amplicons that minimise the longest uncovered stretch.

    ``fixed`` are positions already covered (the marker gene). The fragment is split by
    them into gaps, and each new amplicon goes into whichever gap currently has the
    largest sub-gap (gap / (points already in it + 1)), which minimises the maximum
    sub-gap; within a gap the points are evenly spaced. Without ``fixed`` this is plain
    even spacing (L/4, L/2, 3L/4 for n = 3).
    """
    pts = sorted({0, frag_len, *[int(x) for x in (fixed or []) if 0 < x < frag_len]})
    gaps = list(zip(pts, pts[1:]))
    k = [0] * len(gaps)
    for _ in range(n):
        i = max(range(len(gaps)), key=lambda j: ((gaps[j][1] - gaps[j][0]) / (k[j] + 1), -j))
        k[i] += 1
    centres = []
    for (a, b), kk in zip(gaps, k):
        for j in range(1, kk + 1):
            c = a + (b - a) * j / (kk + 1)
            centres.append(int(min(max(c, margin), frag_len - margin)))
    return sorted(centres)


def _coverage(design: ScreeningDesign, JL: int, frag_len: int) -> None:
    """Record the stretches of the BGC between neighbouring interior amplicons."""
    mids = sorted(int((a.start + a.end) / 2) - JL for a in design.amplicons
                  if a.role in ("marker-gene", "integrity"))
    pts = [0] + [m for m in mids if 0 < m < frag_len] + [frag_len]
    design.coverage_gaps = list(zip(pts, pts[1:]))
    design.max_gap_span = max(design.coverage_gaps, key=lambda g: g[1] - g[0])
    design.max_gap_bp = design.max_gap_span[1] - design.max_gap_span[0]


def design_screening(construct: str, left_junction: int, fragment_len: int,
                     fragment_seq: Optional[str] = None, bgc_record=None,
                     config: Optional[PipelineConfig] = None,
                     host_seq: Optional[str] = None) -> ScreeningDesign:
    """Design the colony-PCR panel for an assembled construct.

    ``left_junction`` is the construct index where the BGC fragment starts (the
    vector/fragment boundary); the right junction is ``left_junction +
    fragment_len``. ``fragment_seq`` and ``bgc_record`` (a GenBank SeqRecord)
    enable marker-gene detection; without them evenly spaced integrity
    amplicons are used instead. ``host_seq`` (e.g. the yeast genome, one string)
    adds a host background check.
    """
    config = config or PipelineConfig()
    construct = _clean(construct)
    L = len(construct)
    JL, JR = left_junction, left_junction + fragment_len
    design = ScreeningDesign(construct_length=L, left_junction=JL, right_junction=JR)
    if not construct or fragment_len <= 0 or JR > L:
        design.notes.append("screening primers not designed: construct/junction coordinates are inconsistent")
        return design
    cons = Background(construct, circular=True)
    host = Background(host_seq, circular=False) if host_seq else None
    if host is None:
        design.notes.append("host genome not provided: primer specificity was checked on the "
                            "construct only (pass --yeast-genome for a host background check)")

    chosen: list = []
    counter = {"i": 0}

    def next_range() -> tuple:
        return size_slot(counter["i"])

    def run(name, role, target, win, lreg, rreg):
        size_range = next_range()
        counter["i"] += 1
        amp = _design_one(name, role, target, construct, cons, win, lreg, rreg,
                          size_range, config, chosen, host)
        if amp is None:
            design.failed.append(f"{name} ({target}): no primer pair satisfied the constraints")
        else:
            design.amplicons.append(amp)
            chosen.append(amp)
        return amp

    # --- junctions: one primer on each side of the vector/BGC boundary ---------
    smax = next_range()[1]
    run("JL", "junction-left", f"vector | BGC left end (construct position {JL + 1})",
        (JL - smax, JL + smax), (JL - smax + 30, JL - 20), (JL + 20, JL + smax - 30))
    smax = next_range()[1]
    run("JR", "junction-right", f"BGC right end | vector (construct position {JR + 1})",
        (JR - smax, JR + smax), (JR - smax + 30, JR - 20), (JR + 20, JR + smax - 30))

    # --- marker gene -----------------------------------------------------------
    has_features = bool(bgc_record is not None and getattr(bgc_record, "features", None))
    if config.screening_marker_gene:
        if not fragment_seq:
            raise ValueError("--marker-gene needs the fragment sequence")
        marker = find_marker_gene(bgc_record, fragment_seq, [], locus_tag=config.screening_marker_gene)
    else:
        marker = (find_marker_gene(bgc_record, fragment_seq, config.screening_marker_keywords)
                  if fragment_seq and has_features else None)
    n_spaced = max(0, config.screening_n_spaced)
    marker_ok, marker_mid = False, None
    if marker:
        design.marker = marker
        gs, ge = JL + marker["start"], JL + marker["end"]
        if marker.get("n_matches", 1) > 1:
            alts = ", ".join(
                f"{a['name']} ({a['start'] / 1000:.1f}-{a['end'] / 1000:.1f} kb"
                + (f", {a['product']}" if a["product"] else "") + ")"
                for a in marker.get("alternatives", [])[:5])
            design.notes.append(
                f"keyword '{marker['keyword']}' matched {marker['n_matches']} CDS; "
                f"{marker['name']} (nearest the BGC centre) was used. Other candidates: {alts}. "
                "Pick one explicitly with --marker-gene LOCUS_TAG")
        how = ("selected with --marker-gene" if marker["keyword"] == "--marker-gene"
               else f"matched '{marker['keyword']}' in: \"{marker.get('evidence', '')}\"")
        amp = run("MK", "marker-gene",
                  f"{marker['name']} ({marker['product'] or marker['keyword']}) at construct "
                  f"{gs + 1}-{ge}; {how}",
                  (gs - 150, ge + 150), (gs, ge - 60), (gs + 60, ge + 60))
        marker_ok = amp is not None
        if marker_ok:
            marker_mid = int((amp.start + amp.end) / 2) - JL
        else:
            design.notes.append(f"marker gene {marker['name']} is too short or has no suitable primers "
                                "for a marker amplicon; a spaced amplicon was added instead")
    elif has_features:
        design.notes.append("no CDS matched the marker keywords "
                            f"{config.screening_marker_keywords}; spaced amplicons used instead "
                            "(use --marker-keyword or --marker-gene to choose a gene)")
    else:
        design.notes.append("BGC has no annotation (FASTA input): no marker gene; spaced integrity "
                            "amplicons used instead (give --genbank to target a core gene)")
    if not marker_ok:
        n_spaced += 1                      # keep the panel size when there is no marker

    # --- integrity amplicons: minimise the longest stretch without an amplicon --
    centres = _allocate_centres(fragment_len, n_spaced, [marker_mid] if marker_ok else [])
    for i, c in enumerate(centres, 1):
        size_range = next_range()
        counter["i"] += 1
        half = size_range[1] // 2 + 150
        amp, used = None, c
        for off in (0, 800, -800, 1600, -1600):          # nudge along the BGC if the spot has no primers
            used = min(max(c + off, 700), fragment_len - 700)
            centre = JL + used
            amp = _design_one(f"IN{i}", "integrity",
                              f"BGC interior, ~{used / 1000:.1f} kb from the left end "
                              f"(construct ~{centre + 1})", construct, cons,
                              (centre - half, centre + half), (centre - half, centre),
                              (centre, centre + half), size_range, config, chosen, host)
            if amp is not None:
                break
        if amp is None:
            design.failed.append(f"IN{i} (BGC interior, ~{c / 1000:.1f} kb): no primer pair "
                                 "satisfied the constraints")
        else:
            if used != c:
                design.notes.append(f"IN{i} was moved {abs(used - c)} bp along the BGC "
                                    "because no primer pair fit at the planned position")
            design.amplicons.append(amp)
            chosen.append(amp)

    # --- coverage of the BGC by interior amplicons -----------------------------
    if design.amplicons:
        _coverage(design, JL, fragment_len)
        if design.max_gap_bp > config.screening_max_gap_bp:
            a0, b0 = design.max_gap_span
            design.notes.append(
                f"largest stretch of the BGC without an interior amplicon: "
                f"{design.max_gap_bp / 1000:.1f} kb ({a0 / 1000:.1f}-{b0 / 1000:.1f} kb from the "
                "left end); raise --screening-spaced for denser coverage")

    # --- panel-level checks: multiplex cross products ---------------------------
    if len(design.amplicons) > 1:
        all_primers = {}
        for a in design.amplicons:
            all_primers[a.forward.name] = a.forward.sequence
            all_primers[a.reverse.name] = a.reverse.sequence
        own = {(a.forward.name, a.reverse.name) for a in design.amplicons}
        # Any forward primer pairs with any reverse primer further downstream, but only
        # products close to the intended sizes compete or confuse the gel; long ones don't.
        limit = max(a.product_size for a in design.amplicons) + 500
        for p in predict_amplicons(all_primers, construct, circular=True, max_size=limit):
            if (p["forward"], p["reverse"]) not in own:
                design.multiplex_products.append(p)
        if design.multiplex_products:
            design.notes.append(
                f"{len(design.multiplex_products)} cross-pair product(s) up to {limit} bp predicted "
                "if all primers are combined in one tube (sizes: "
                + ", ".join(str(p['size']) for p in design.multiplex_products[:6])
                + " bp); run the pairs as separate reactions or check the band sizes")
    return design


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------
def screening_rows(design: ScreeningDesign) -> list:
    rows = []
    for a in design.amplicons:
        for p, partner in ((a.forward, a.reverse), (a.reverse, a.forward)):
            rows.append({
                "primer": p.name, "amplicon": a.name, "role": a.role, "target": a.target,
                "sequence_5to3": p.sequence, "strand": p.strand, "length": p.length,
                "construct_start_1based": p.start + 1, "construct_end": p.end,
                "tm_c": p.tm, "gc_percent": p.gc_percent, "partner": partner.name,
                "product_bp": a.product_size, "design_level": a.level,
                "issues": "; ".join(a.issues), "warnings": "; ".join(a.warnings),
            })
    return rows


def export_screening(design: ScreeningDesign, path: str) -> str:
    """Write the screening primers to a CSV ready for a synthesis order."""
    rows = screening_rows(design)
    cols = ["primer", "amplicon", "role", "target", "sequence_5to3", "strand", "length",
            "construct_start_1based", "construct_end", "tm_c", "gc_percent", "partner",
            "product_bp", "design_level", "issues", "warnings"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return path


__all__ = [
    "SIZE_LADDER", "MIN_PRODUCT_SEPARATION", "size_slot", "find_primer_sites", "predict_amplicons", "find_marker_gene",
    "design_screening", "screening_rows", "export_screening",
]
