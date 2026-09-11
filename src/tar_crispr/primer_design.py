"""Step 5: Tailed primer design for the capture vector.

A tailed primer consists of:
    [homology arm (40-60 nt)] + [vector annealing region (18-25 nt)]

The homology arm enables recombination with the BGC fragment, while the
annealing region provides PCR-specific priming on the linearized vector.
primer3-py is used to optimize the annealing region for Tm, hairpins,
dimers, and GC clamp.
"""
import csv
import warnings
from typing import Optional

import primer3

from tar_crispr.config import PipelineConfig, TailedPrimer, HomologyArm


warnings.filterwarnings("ignore", category=UserWarning)


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTNacgtn", "TGCANTGCAN")
    return seq.translate(comp)[::-1]


def design_tailed_primer(arm: HomologyArm,
                         vector_seq: str,
                         vector_cut_position: int,
                         end: str,
                         config: PipelineConfig) -> TailedPrimer:
    """Design a tailed primer for one end of the capture vector.

    Parameters
    ----------
    arm : HomologyArm
        The validated homology arm for this end.
    vector_seq : str
        The full capture vector sequence.
    vector_cut_position : int
        0-based position where the vector is linearized (for this end).
    end : str
        ``"left"`` or ``"right"``.
    config : PipelineConfig
        Pipeline configuration.

    Returns
    -------
    TailedPrimer
        The designed primer with all metrics.
    """
    vec = vector_seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")

    if end == "left":
        anneal_end = vector_cut_position
        annealing_region = _optimize_annealing_length(
            vec, anneal_end,
            config.primer_anneal_min, config.primer_anneal_max,
            config.primer_tm_min, config.primer_tm_max,
            max_shift=100,
        )
        primer_type = "FORWARD"
    else:
        rc_vec = _reverse_complement(vec)
        rc_cut = len(vec) - vector_cut_position
        anneal_end = rc_cut
        annealing_region = _optimize_annealing_length(
            rc_vec, anneal_end,
            config.primer_anneal_min, config.primer_anneal_max,
            config.primer_tm_min, config.primer_tm_max,
            max_shift=100,
        )
        primer_type = "REVERSE"

    full_primer = arm.sequence + annealing_region
    # Tm is calculated for the annealing region (the part that binds the template
    # during PCR extension). The tail provides homology but does not anneal.
    tm = _calculate_tm(annealing_region) if annealing_region else 0.0
    gc = _gc_percent(annealing_region)
    # Hairpin: check the 3' end of the annealing region (last 20 bp max)
    hp_check_seq = annealing_region[-20:] if len(annealing_region) > 20 else annealing_region
    hp = _check_hairpin(hp_check_seq)
    # Self-dimer: this primer's annealing region against itself
    # (informational — real PCR failure usually comes from cross-dimer)
    self_dimer = _check_dimer(annealing_region, annealing_region)
    # Palindromic 3' end is a softer signal than full self-dimer
    palindromic_3p = _has_palindromic_3p(annealing_region)

    # issues   = hard failures (block synthesis / PCR)
    # warnings = soft signals (informational, primer is still usable)
    issues = []
    warnings = []
    if len(annealing_region) == 0:
        issues.append("Annealing region is empty — check vector cut positions")
    else:
        if not _check_gc_clamp(annealing_region):
            warnings.append("No GC clamp at 3' end of annealing region")
        if tm < config.primer_tm_min:
            # Relaxed tolerance: accept down to tm_min - 3 C as a warning
            if tm < config.primer_tm_min - 3.0:
                issues.append(f"Tm {tm:.1f}C too far below minimum {config.primer_tm_min}C")
            else:
                warnings.append(f"Tm {tm:.1f}C below target minimum {config.primer_tm_min}C (within relaxed tolerance)")
        elif tm > config.primer_tm_max:
            if tm > config.primer_tm_max + 3.0:
                issues.append(f"Tm {tm:.1f}C too far above maximum {config.primer_tm_max}C")
            else:
                warnings.append(f"Tm {tm:.1f}C above target maximum {config.primer_tm_max}C (within relaxed tolerance)")
    if hp:
        warnings.append("Significant hairpin predicted (self)")
    if self_dimer:
        warnings.append("Significant self-dimer predicted")
    if palindromic_3p:
        warnings.append("Palindromic 3' end of annealing region")

    valid = len(issues) == 0

    return TailedPrimer(
        name=f"{primer_type}_{end}",
        sequence=full_primer,
        tail=arm.sequence,
        annealing_region=annealing_region,
        tm=round(tm, 1),
        gc_percent=round(gc, 1),
        hairpin=hp,
        self_dimer=self_dimer,
        cross_dimer=False,   # set later in design_tailed_primers
        issues=issues,
        warnings=warnings,
        valid=valid,
    )


def _has_palindromic_3p(seq: str, min_len: int = 4) -> bool:
    """Return True if the 3' end of *seq* contains a palindromic motif.

    Common PCR-unfriendly 3' endings: CGCG, GCGC, ATAT, TATA, AATT, TTAA,
    GATC, CATG, etc. We check for the canonical 4 bp palindromes.
    """
    if len(seq) < min_len:
        return False
    tail = seq[-min_len:].upper()
    palindromes = {"CGCG", "GCGC", "ATAT", "TATA", "AATT", "TTAA",
                   "GATC", "CATG", "ACGT", "TGCA"}
    return tail in palindromes


def _optimize_annealing_length(vector_seq: str,
                               search_end: int,
                               min_len: int,
                               max_len: int,
                               tm_min: float,
                               tm_max: float,
                               max_shift: int = 20) -> str:
    """Find the annealing region whose Tm falls within [tm_min, tm_max].

    Slides a window along the vector within +/- max_shift bp of the cut site,
    trying all lengths from max_len down to min_len. Prefers candidates that
    are in-range AND have a GC clamp; otherwise falls back to the candidate
    with Tm closest to the interval midpoint.
    """
    target_tm = (tm_min + tm_max) / 2
    # Relaxed acceptance window: some vectors (e.g. pCAP03 with SwaI) have
    # cut sites in AT-rich terminators, where strict Tm 58-62 C is unattainable.
    # We accept 55-65 C as "good enough" and 50-55 C as "provisional".
    relaxed_min = tm_min - 3.0
    relaxed_max = tm_max + 3.0

    best_seq = None
    best_score = float("inf")
    best_info = ("", -1, -1)  # (kind, shift, length)

    # IMPORTANT: only allow sliding UPSTREAM (negative d) from search_end.
    # A positive d would place the annealing past the cut site, which is
    # biologically wrong (the primer would not be adjacent to the cut).
    # We permit at most 3 bp of slack on the downstream side.
    downstream_slack = 3

    for length in range(max_len, min_len - 1, -1):
        for shift in range(0, max_shift + 1):
            if shift == 0:
                shifts = (0,)
            else:
                shifts = (-shift,)
                if shift <= downstream_slack:
                    shifts = (shift,) + shifts  # prefer upstream first
            for d in shifts:
                end = search_end + d
                start = end - length
                if start < 0 or end > len(vector_seq):
                    continue
                candidate = vector_seq[start:end]
                if len(candidate) != length:
                    continue
                tm = _calculate_tm(candidate)
                has_clamp = _check_gc_clamp(candidate)

                if tm_min <= tm <= tm_max and has_clamp:
                    return candidate  # strict hit, closest to cut wins

                # remember relaxed hits, preferring those closer to cut
                if relaxed_min <= tm <= relaxed_max and has_clamp:
                    if best_info[0] != "strict" and (
                        best_info[0] == "" or shift < best_info[1]
                    ):
                        best_seq = candidate
                        best_info = ("relaxed", shift, length)
                    continue

                # Fallback scoring: prefer close-to-target AND close to cut
                score = abs(tm - target_tm)
                if not (relaxed_min <= tm <= relaxed_max):
                    score += 100
                score += shift * 1.0  # strong penalty for distance
                if score < best_score:
                    best_score = score
                    best_seq = candidate
                    best_info = ("fallback", shift, length)

    if best_seq is None:
        best_seq = vector_seq[max(0, search_end - 20):search_end]
    return best_seq



def _calculate_tm(seq: str) -> float:
    """Calculate Tm (Celsius) using primer3's nearest-neighbor model."""
    if not seq or len(seq) == 0:
        return 0.0
    try:
        return primer3.calc_tm(seq)
    except Exception:
        return _wallace_tm(seq)


def _wallace_tm(seq: str) -> float:
    """Wallace rule Tm estimation (fallback)."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in "GC")
    if len(seq) < 14:
        return float(len(seq) * 2 + gc * 2)
    return round(64.9 + 41 * (gc - 16.4) / len(seq), 1)


def _gc_percent(seq: str) -> float:
    """Return GC percentage (0-100)."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return gc / len(seq) * 100


def _check_hairpin(seq: str) -> bool:
    """Check for significant hairpin formation (ΔG < -2.0 kcal/mol at 37 C).

    primer3 limits thermodynamic analysis to 60 bp, so we analyze the 3'
    half of the primer (most relevant for hairpin-induced PCR artifacts).
    """
    if len(seq) > 60:
        seq = seq[-60:]
    try:
        result = primer3.calc_hairpin(seq, temp_c=37)
        return result.dg < -2.0
    except Exception:
        return _check_hairpin_heuristic(seq)


def _check_hairpin_heuristic(seq: str) -> bool:
    """Fallback: detect 4+ bp inverted repeats that could form a hairpin."""
    n = len(seq)
    for stem_len in range(4, n // 2 + 1):
        for i in range(n - stem_len * 2 - 2):
            seg = seq[i:i + stem_len]
            rc_seg = _reverse_complement(seg)
            for j in range(i + stem_len + 3, n - stem_len + 1):
                if seq[j:j + stem_len] == rc_seg:
                    return True
    return False


def _check_dimer(seq: str, other: str) -> bool:
    """Check for primer-dimer formation between two sequences."""
    # Use the shorter region for dimer analysis (primer3 limit: 60 bp)
    s1 = seq[-60:] if len(seq) > 60 else seq
    s2 = other[-60:] if len(other) > 60 else other
    try:
        result = primer3.calc_heterodimer(s1, s2, temp_c=37)
        return result.dg < -2.0
    except Exception:
        return _check_dimer_heuristic(s1, s2)


def _check_dimer_heuristic(seq: str, other: str) -> bool:
    """Fallback: detect 4+ bp 3' complementarity between two primers."""
    for overlap in range(4, min(len(seq), len(other)) + 1):
        if seq[-overlap:].upper() == _reverse_complement(other[:overlap].upper()):
            return True
        if seq[:overlap].upper() == _reverse_complement(other[-overlap:].upper()):
            return True
    return False


def _check_gc_clamp(seq: str, min_gc: int = 2) -> bool:
    """Check for GC clamp: at least *min_gc* G/C in the last 5 bases."""
    if len(seq) < 5:
        return False
    last_5 = seq[-5:].upper()
    return sum(1 for b in last_5 if b in "GC") >= min_gc


def design_tailed_primers(left_arm: HomologyArm,
                          right_arm: HomologyArm,
                          vector_seq: str,
                          vector_cut_left: int,
                          vector_cut_right: int,
                          config: Optional[PipelineConfig] = None) -> dict:
    """Design tailed primers for both ends of the capture vector.

    Parameters
    ----------
    left_arm : HomologyArm
        Left homology arm from the BGC fragment.
    right_arm : HomologyArm
        Right homology arm from the BGC fragment.
    vector_seq : str
        Capture vector sequence.
    vector_cut_left : int
        0-based left linearization coordinate.
    vector_cut_right : int
        0-based right linearization coordinate.
    config : PipelineConfig, optional
        Pipeline configuration (defaults if omitted).

    Returns
    -------
    dict
        ``{"left": TailedPrimer, "right": TailedPrimer}``
    """
    if config is None:
        config = PipelineConfig()

    left_primer = design_tailed_primer(left_arm, vector_seq,
                                       vector_cut_left, "left", config)
    right_primer = design_tailed_primer(right_arm, vector_seq,
                                        vector_cut_right, "right", config)

    tm_delta = abs(left_primer.tm - right_primer.tm)
    # Tm delta is a hard failure only if it's large; small differences are warnings.
    if tm_delta > config.primer_max_delta_tm + 1.0:
        left_primer.issues.append(
            f"Tm difference {tm_delta:.1f}C exceeds max {config.primer_max_delta_tm}C"
        )
        right_primer.issues.append(
            f"Tm difference {tm_delta:.1f}C exceeds max {config.primer_max_delta_tm}C"
        )
        left_primer.valid = False
        right_primer.valid = False
    elif tm_delta > config.primer_max_delta_tm:
        msg = f"Tm difference {tm_delta:.1f}C slightly exceeds max {config.primer_max_delta_tm}C"
        left_primer.warnings.append(msg)
        right_primer.warnings.append(msg)

    # Cross-dimer between forward and reverse annealing regions.
    # We compute 3'-end complementarity explicitly, because primer3's
    # heterodimer check flags many harmless internal pairings.
    cross_dimer = _has_3p_cross_complementarity(
        left_primer.annealing_region, right_primer.annealing_region
    )
    left_primer.cross_dimer = cross_dimer
    right_primer.cross_dimer = cross_dimer
    if cross_dimer:
        msg = "3' cross-dimer detected between forward and reverse"
        left_primer.warnings.append(msg)
        right_primer.warnings.append(msg)

    return {"left": left_primer, "right": right_primer}


def _has_3p_cross_complementarity(seq1: str, seq2: str,
                                  min_len: int = 4) -> bool:
    """Check whether the 3' end of one primer can anneal to the 3' end of the other.

    This is the biologically meaningful cross-dimer check for PCR: if either
    3' end is complementary to the other 3' end for >= min_len bp, the pair
    can form primer-dimers that get extended by the polymerase.

    Internal pairings (5' with 5', or 5' with 3') are much less harmful and
    are intentionally NOT flagged here.
    """
    if not seq1 or not seq2:
        return False
    s1 = seq1.upper()
    s2 = seq2.upper()
    max_ov = min(len(s1), len(s2), 12)  # cap search at 12 bp overlap

    for ov in range(max_ov, min_len - 1, -1):
        # 3' of s1 vs 3' of s2 (reverse-complement of one must match)
        if _reverse_complement(s1[-ov:]) == s2[-ov:]:
            return True
        # 3' of s1 vs 5' of s2 (asymmetric dimer)
        if _reverse_complement(s1[-ov:]) == s2[:ov]:
            return True
        # 5' of s1 vs 3' of s2
        if _reverse_complement(s1[:ov]) == s2[-ov:]:
            return True
    return False


def export_primers(primers: dict, path: str) -> str:
    """Export primers to a CSV file for synthesis ordering.

    Parameters
    ----------
    primers : dict
        Output from ``design_tailed_primers``.
    path : str
        Output CSV file path.

    Returns
    -------
    str
        Path to the written CSV file.
    """
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "Primer_Name", "Sequence", "Tail_Length", "Annealing_Length",
            "Tm_C", "GC_Percent", "Hairpin", "Self_Dimer", "Cross_Dimer",
            "Valid", "Issues", "Warnings"
        ])
        for key in ("left", "right"):
            p = primers[key]
            writer.writerow([
                p.name, p.sequence, len(p.tail), len(p.annealing_region),
                p.tm, p.gc_percent, p.hairpin, p.self_dimer, p.cross_dimer,
                p.valid,
                "; ".join(p.issues) if p.issues else "",
                "; ".join(p.warnings) if p.warnings else "",
            ])
    return path


__all__ = [
    "design_tailed_primer",
    "design_tailed_primers",
    "export_primers",
]
