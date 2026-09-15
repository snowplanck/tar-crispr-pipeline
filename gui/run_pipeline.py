"""Thin wrapper around the TAR-CRISPR pipeline for the Streamlit GUI.

This module contains no pipeline logic of its own: it just calls the
existing functions in ``src/tar_crispr/`` in the same order as
``cli.py`` does, and returns a dict with everything the GUI needs to
render the report.
"""
from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path
from typing import Callable, Optional

from tar_crispr.config import PipelineConfig
from tar_crispr.sequence_io import (
    read_sequence,
    validate_sequence,
    extract_cluster_bounds,
    get_flanking_sequence,
    check_external_tools,
)
from tar_crispr.pam_finder import design_sgRNAs
from tar_crispr.fragment_ends import extract_fragment
from tar_crispr.homology_arms import design_homology_arms
from tar_crispr.primer_design import design_tailed_primers, export_primers
from tar_crispr.assembly_sim import simulate_pydna_assembly
from tar_crispr.report_generator import generate_report


ProgressCB = Optional[Callable[[str, float], None]]


def run_pipeline(
    bgc_path: str | Path,
    vector_path: str | Path,
    genome_path: str | Path,
    output_dir: str | Path,
    start: Optional[int] = None,
    end: Optional[int] = None,
    vector_enzyme: Optional[str] = None,
    use_blast: bool = False,
    pam_window: int = 500,
    arm_length: int = 50,
    top_n: int = 5,
    tm_min: float = 58.0,
    tm_max: float = 62.0,
    max_mismatches: int = 3,
    progress: ProgressCB = None,
) -> dict:
    """Run the full pipeline and return artifacts for rendering.

    Parameters mirror the CLI flags. ``progress`` is an optional callable
    invoked as ``progress(message, fraction)`` after each step, with
    ``fraction`` in [0, 1].
    """
    def _step(msg: str, frac: float) -> None:
        if progress is not None:
            progress(msg, frac)

    bgc_path = Path(bgc_path)
    vector_path = Path(vector_path)
    genome_path = Path(genome_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Config ---
    config = PipelineConfig(
        pam_window=pam_window,
        homology_arm_length=arm_length,
        top_n_sgRNAs=top_n,
        primer_tm_min=tm_min,
        primer_tm_max=tm_max,
        max_mismatches=max_mismatches,
    )
    config = check_external_tools(config)
    if not use_blast:
        config = replace(config, blast_available=False)

    # --- Step 1: read sequences ---
    _step("Reading input sequences...", 0.05)
    bgc_record = read_sequence(bgc_path)
    vector_record = read_sequence(vector_path)
    genome_record = read_sequence(genome_path)

    bgc_stats = validate_sequence(bgc_record)
    vector_stats = validate_sequence(vector_record)

    if not bgc_stats["valid"]:
        raise ValueError(
            f"BGC has invalid characters: {bgc_stats['invalid_chars']}"
        )
    if not vector_stats["valid"]:
        raise ValueError(
            f"Vector has invalid characters: {vector_stats['invalid_chars']}"
        )

    # --- Step 2: cluster coordinates ---
    _step("Determining cluster bounds...", 0.15)
    has_features = bool(getattr(bgc_record, "features", None))
    if (start is None or end is None) and not has_features:
        raise ValueError(
            "Cannot determine cluster bounds: the BGC file has no annotated "
            "features and --start/--end were not provided."
        )
    cluster = extract_cluster_bounds(
        genome_record if has_features is False else bgc_record,
        start, end,
    ) if (start is not None and end is not None) else extract_cluster_bounds(
        genome_record, start, end
    )

    # --- Step 3: flanks ---
    _step("Extracting flanks...", 0.25)
    upstream, cluster_seq, downstream = get_flanking_sequence(
        genome_record, cluster, config.pam_window
    )

    # --- Step 4: sgRNAs ---
    _step("Designing sgRNAs...", 0.40)
    genome_seq = str(genome_record.seq).upper()
    sgRNAs = design_sgRNAs(cluster_seq, upstream, downstream, genome_seq, config)
    if not sgRNAs["left"] or not sgRNAs["right"]:
        raise ValueError(
            "No sgRNA candidates found at one or both boundaries. "
            "Try increasing the PAM window or check the coordinates."
        )

    selected = {"left": sgRNAs["left"][0], "right": sgRNAs["right"][0]}

    # --- Step 5: fragment ---
    _step("Excising fragment...", 0.55)
    fragment = extract_fragment(genome_seq, selected["left"], selected["right"])

    # --- Step 6: homology arms ---
    _step("Designing homology arms...", 0.65)
    arms = design_homology_arms(
        fragment.sequence, genome_seq, str(vector_record.seq).upper(), config
    )

    # --- Step 7: primers ---
    _step("Designing tailed primers...", 0.75)
    vector_seq = str(vector_record.seq).upper()
    cut_left, cut_right = _resolve_vector_cuts(
        vector_seq, vector_enzyme, config
    )
    primers = design_tailed_primers(
        arms["left"], arms["right"], vector_seq,
        cut_left, cut_right, config,
    )

    # --- Step 8: assembly ---
    _step("Simulating assembly...", 0.85)
    assembly = simulate_pydna_assembly(
        fragment.sequence, vector_seq,
        arms["left"], arms["right"], cut_left, cut_right,
    )

    # --- Step 9: report ---
    _step("Generating report...", 0.95)
    generate_report(
        sgRNAs, fragment.sequence, arms, primers, assembly,
        cluster, config, bgc_stats, str(output_dir),
    )
    # report_generator writes report.md / report.html / cluster_map.svg
    # and cli.py writes primers.csv via export_primers. Do it here too.
    export_primers(primers, str(output_dir / "primers.csv"))

    _step("Done.", 1.0)

    return {
        "sgRNAs": sgRNAs,
        "selected": selected,
        "fragment": fragment,
        "arms": arms,
        "primers": primers,
        "assembly": assembly,
        "cluster": cluster,
        "config": config,
        "bgc_stats": bgc_stats,
        "vector_stats": vector_stats,
        "vector_cuts": (cut_left, cut_right),
        "output_dir": output_dir,
        "report_md": output_dir / "report.md",
        "report_html": output_dir / "report.html",
        "cluster_map_svg": output_dir / "cluster_map.svg",
        "primers_csv": output_dir / "primers.csv",
    }


def _resolve_vector_cuts(
    vector_seq: str,
    enzyme: Optional[str],
    config: PipelineConfig,
) -> tuple[int, int]:
    """Return (cut_left, cut_right) 0-based.

    If ``enzyme`` is given, locate its single recognition site; otherwise
    fall back to defaults stored in the config.
    """
    if enzyme:
        site = _enzyme_site(enzyme)
        if site is None:
            raise ValueError(f"Unknown enzyme: {enzyme}")
        positions = _find_all(vector_seq, site)
        if len(positions) != 1:
            raise ValueError(
                f"{enzyme} has {len(positions)} sites in the vector; "
                f"expected exactly one."
            )
        pos = positions[0]  # 0-based start of site
        # cut in the middle of the site
        mid = pos + len(site) // 2
        return (mid - 4, mid + 4)
    return (0, 8)  # fallback: user must override via cut-left/right


def _find_all(haystack: str, needle: str) -> list[int]:
    hits, start = [], 0
    while True:
        i = haystack.find(needle, start)
        if i == -1:
            return hits
        hits.append(i)
        start = i + 1


def _enzyme_site(name: str) -> Optional[str]:
    table = {
        "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
        "XhoI": "CTCGAG", "SalI": "GTCGAC", "NotI": "GCGGCCGC",
        "SwaI": "ATTTAAAT", "PmeI": "GTTTAAAC", "NruI": "TCGCGA",
        "AatII": "GACGTC", "PvuI": "CGATCG", "KpnI": "GGTACC",
        "NheI": "GCTAGC", "ClaI": "ATCGAT", "BglII": "AGATCT",
        "ScaI": "AGTACT", "SpeI": "ACTAGT",
    }
    return table.get(name)

