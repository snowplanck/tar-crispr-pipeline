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
        anneal_start = max(0, vector_cut_position - config.primer_anneal_max)
        anneal_end = vector_cut_position
        annealing_region = _optimize_annealing_length(
            vec, anneal_start, anneal_end,
            config.primer_anneal_min, config.primer_anneal_max,
            config.primer_tm_min, config.primer_tm_max,
        )
        primer_type = "FORWARD"
    else:
        rc_vec = _reverse_complement(vec)
        rc_cut = len(vec) - vector_cut_position
        anneal_start = max(0, rc_cut - config.primer_anneal_max)
        anneal_end = rc_cut
        annealing_region = _optimize_annealing_length(
            rc_vec, anneal_start, anneal_end,
            config.primer_anneal_min, config.primer_anneal_max,
            config.primer_tm_min, config.primer_tm_max,
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
    # Self-dimer: check the annealing region only
    dimer = _check_dimer(annealing_region, annealing_region)

    issues = []
    if not _check_gc_clamp(annealing_region):
        issues.append("No GC clamp at 3' end of annealing region")
    if len(annealing_region) == 0:
        issues.append("Annealing region is empty — check vector cut positions")
    elif tm < config.primer_tm_min:
        issues.append(f"Tm {tm:.1f}C below minimum {config.primer_tm_min}C")
    elif tm > config.primer_tm_max:
        issues.append(f"Tm {tm:.1f}C above maximum {config.primer_tm_max}C")
    if hp:
        issues.append("Significant hairpin predicted")
    if dimer:
        issues.append("Significant dimer predicted")

    valid = len(issues) == 0

    return TailedPrimer(
        name=f"{primer_type}_{end}",
        sequence=full_primer,
        tail=arm.sequence,
        annealing_region=annealing_region,
        tm=round(tm, 1),
        gc_percent=round(gc, 1),
        hairpin=hp,
        dimer=dimer,
        issues=issues,
        valid=valid,
    )


def _optimize_annealing_length(vector_seq: str,
                               search_start: int,
                               search_end: int,
                               min_len: int,
                               max_len: int,
                               tm_min: float,
                               tm_max: float) -> str:
    """Find the annealing region whose Tm falls within [tm_min, tm_max].

    Tries lengths from *max_len* down to *min_len*. Returns the first candidate
    whose Tm is in range, otherwise the one closest to the midpoint.
    """
    target_tm = (tm_min + tm_max) / 2
    best_seq = None
    best_delta = float("inf")

    for length in range(max_len, min_len - 1, -1):
        s = max(0, search_end - length)
        candidate = vector_seq[s:search_end]
        if len(candidate) != length:
            continue
        tm = _calculate_tm(candidate)
        delta = abs(tm - target_tm)
        if tm_min <= tm <= tm_max:
            return candidate
        if delta < best_delta:
            best_delta = delta
            best_seq = candidate

    if best_seq is None:
        s = max(0, search_end - 20)
        best_seq = vector_seq[s:search_end]
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
    if tm_delta > config.primer_max_delta_tm:
        left_primer.issues.append(
            f"Tm difference {tm_delta:.1f}C exceeds max {config.primer_max_delta_tm}C"
        )
        left_primer.valid = False

    # Check cross-dimer between forward and reverse annealing regions
    cross_dimer = _check_dimer(left_primer.annealing_region,
                               right_primer.annealing_region)
    if cross_dimer:
        left_primer.issues.append("Cross-dimer detected between forward and reverse")
        left_primer.valid = False

    return {"left": left_primer, "right": right_primer}


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
            "Tm_C", "GC_Percent", "Hairpin", "Dimer", "Valid", "Notes"
        ])
        for key in ("left", "right"):
            p = primers[key]
            writer.writerow([
                p.name, p.sequence, len(p.tail), len(p.annealing_region),
                p.tm, p.gc_percent, p.hairpin, p.dimer, p.valid,
                "; ".join(p.issues) if p.issues else "",
            ])
    return path


__all__ = [
    "design_tailed_primer",
    "design_tailed_primers",
    "export_primers",
]
