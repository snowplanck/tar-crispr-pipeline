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
  partial captures that still have both junctions.

Every amplicon gets a different product size so the panel can be run as one
multiplex or read from a single gel. Candidate pairs come from primer3; each is
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

# Product-size ladder: ranges are >= 150 bp apart so bands resolve on a gel.
SIZE_LADDER = [(300, 450), (500, 650), (700, 850), (900, 1050),
               (1100, 1250), (1300, 1450), (1500, 1650), (1700, 1850)]

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
    parts = []
    for key in ("product", "gene", "gene_functions", "sec_met_domain", "note", "function"):
        parts.extend(str(v) for v in feature.qualifiers.get(key, []))
    return " ".join(parts).lower()


def find_marker_gene(bgc_record, fragment_seq: str, keywords: list) -> Optional[dict]:
    """Locate a core biosynthetic gene of the BGC annotation inside the fragment.

    Keywords are tried in priority order against each CDS's product/gene/
    gene_functions/sec_met_domain/note text. The gene is located by searching its
    nucleotide sequence in the fragment (either strand), so no coordinate
    translation between the annotation record, the genome and the fragment is
    needed. Among CDSs matching the first productive keyword, the one nearest
    the fragment centre is used. Returns ``None`` if nothing matches.
    """
    if bgc_record is None or not getattr(bgc_record, "features", None):
        return None
    fragment = _clean(fragment_seq)
    centre = len(fragment) / 2
    cds = [f for f in bgc_record.features if f.type == "CDS"]
    for kw in keywords:
        kw = kw.lower()
        hits = []
        for f in cds:
            if kw not in _feature_text(f):
                continue
            try:
                gene = _clean(str(f.extract(bgc_record.seq)))
            except Exception:
                continue
            if len(gene) < 90:
                continue
            pos = fragment.find(gene)
            strand = "+"
            if pos < 0:
                pos = fragment.find(_rc(gene))
                strand = "-"
            if pos < 0 or fragment.count(gene if strand == "+" else _rc(gene)) != 1:
                continue
            name = (f.qualifiers.get("locus_tag") or f.qualifiers.get("gene")
                    or f.qualifiers.get("protein_id") or ["CDS"])[0]
            product = (f.qualifiers.get("product") or [""])[0]
            hits.append({"name": str(name), "product": str(product), "keyword": kw,
                         "start": pos, "end": pos + len(gene), "strand": strand})
        if hits:
            return min(hits, key=lambda h: abs((h["start"] + h["end"]) / 2 - centre))
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
def _spaced_centres(frag_len: int, n: int, avoid: Optional[tuple]) -> list:
    centres = []
    for i in range(1, n + 1):
        c = int(frag_len * i / (n + 1))
        if avoid and abs(c - (avoid[0] + avoid[1]) / 2) < 1500:
            shift = 1500 if c >= (avoid[0] + avoid[1]) / 2 else -1500
            c = min(max(c + shift, 700), frag_len - 700)
        centres.append(c)
    return centres


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
        return SIZE_LADDER[min(counter["i"], len(SIZE_LADDER) - 1)]

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
    marker = (find_marker_gene(bgc_record, fragment_seq, config.screening_marker_keywords)
              if fragment_seq and has_features else None)
    n_spaced = max(0, config.screening_n_spaced)
    marker_ok = False
    if marker:
        design.marker = marker
        gs, ge = JL + marker["start"], JL + marker["end"]
        amp = run("MK", "marker-gene",
                  f"{marker['name']} ({marker['product'] or marker['keyword']}) at construct "
                  f"{gs + 1}-{ge}",
                  (gs - 150, ge + 150), (gs, ge - 60), (gs + 60, ge + 60))
        marker_ok = amp is not None
        if not marker_ok:
            design.notes.append(f"marker gene {marker['name']} is too short or has no suitable primers "
                                "for a marker amplicon; a spaced amplicon was added instead")
    elif has_features:
        design.notes.append("no CDS matched the marker keywords "
                            f"{config.screening_marker_keywords}; spaced amplicons used instead "
                            "(use --marker-keyword to choose a gene)")
    else:
        design.notes.append("BGC has no annotation (FASTA input): no marker gene; spaced integrity "
                            "amplicons used instead (give --genbank to target a core gene)")
    if not marker_ok:
        n_spaced += 1                      # keep the panel size when there is no marker

    # --- integrity amplicons spread over the fragment --------------------------
    avoid = (marker["start"], marker["end"]) if marker else None
    for i, c in enumerate(_spaced_centres(fragment_len, n_spaced, avoid), 1):
        size_range = next_range()
        centre = JL + c
        half = size_range[1] // 2 + 150
        run(f"IN{i}", "integrity",
            f"BGC interior, ~{c / 1000:.1f} kb from the left end (construct ~{centre + 1})",
            (centre - half, centre + half), (centre - half, centre), (centre, centre + half))

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
    "SIZE_LADDER", "find_primer_sites", "predict_amplicons", "find_marker_gene",
    "design_screening", "screening_rows", "export_screening",
]
