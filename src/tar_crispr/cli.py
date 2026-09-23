"""CLI entry point for the TAR-CRISPR pipeline.

Usage:
    pipeline-tar-crispr --bgc bgc.fasta --vector vector.fasta --output out/
    pipeline-tar-crispr --bgc bgc.fasta --vector vector.fasta --genome genome.fasta --genbank bgc.gbk --output out/
"""
import sys
from pathlib import Path
from typing import Optional, List

import typer

from tar_crispr.config import PipelineConfig, ClusterInfo, PAMCandidate, HomologyArm
from tar_crispr.sequence_io import (
    read_sequence, read_genome, locate_bgc_in_genome, translate_to_flat,
    validate_sequence, extract_cluster_bounds,
    get_flanking_sequence, check_external_tools,
)
from tar_crispr.pam_finder import design_sgRNAs, find_pam_sites, rank_sgRNAs
from tar_crispr.fragment_ends import extract_fragment
from tar_crispr.homology_arms import design_homology_arms
from tar_crispr.primer_design import design_tailed_primers
from tar_crispr.assembly_sim import simulate_pydna_assembly
from tar_crispr.report_generator import generate_report, generate_construct_map

app = typer.Typer(
    help="TAR-CRISPR: Automated sgRNA + homology arm + primer design for CATCH-style BGC capture.",
    add_completion=False,
)


@app.command("run")
def run(
    bgc: str = typer.Option(..., "--bgc", "-b", help="BGC FASTA/GenBank file"),
    vector: str = typer.Option(..., "--vector", "-v", help="Capture vector FASTA file"),
    output: str = typer.Option("output/", "--output", "-o", help="Output directory"),
    genome: Optional[str] = typer.Option(None, "--genome", "-g", help="Full genome FASTA for specificity"),
    genbank: Optional[str] = typer.Option(None, "--genbank", help="GenBank file with annotations"),
    start: Optional[int] = typer.Option(None, "--start", help="BGC start coordinate (0-based)"),
    end: Optional[int] = typer.Option(None, "--end", help="BGC end coordinate (0-based)"),
    gene_kinds: Optional[List[str]] = typer.Option(
        None, "--gene-kinds",
        help=(
            "Filter antiSMASH region bounds by gene_kind "
            "(e.g. biosynthetic,biosynthetic-additional). "
            "Pass 'all' to use the region span as-is."
        ),
    ),
    vector_cut_left: int = typer.Option(None, "--vector-cut-left", help="Vector left linearization coordinate (0-based)"),
    vector_cut_right: Optional[int] = typer.Option(None, "--vector-cut-right", help="Vector right linearization coordinate (0-based)"),
    vector_enzyme: str = typer.Option("EcoRI", "--vector-enzyme", help="Restriction enzyme for vector linearization (when cut coords not specified)"),
    top_n: int = typer.Option(5, "--top-n", help="Top N sgRNAs per end"),
    arm_length: int = typer.Option(50, "--arm-length", help="Homology arm length (bp)"),
    pam_window: int = typer.Option(500, "--pam-window", help="Flanking window for PAM search (bp)"),
    max_shift: int = typer.Option(100, "--max-shift", help="Max shift for arm optimization (bp)"),
    max_mismatches: int = typer.Option(3, "--max-mismatches", help="Max mismatches for specificity"),
    tm_min: float = typer.Option(58.0, "--tm-min", help="Minimum primer Tm"),
    tm_max: float = typer.Option(62.0, "--tm-max", help="Maximum primer Tm"),
    avoid_enzymes: Optional[List[str]] = typer.Option(None, "--avoid-enzyme", help="Restriction enzymes to avoid in arms"),
    force_fallback: bool = typer.Option(False, "--no-blast", help="Force k-mer fallback instead of BLAST"),
    auto_select: bool = typer.Option(True, "--auto-select/--manual-select", help="Auto-select top-ranked sgRNA per end"),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Verbose output"),
):
    """Run the full TAR-CRISPR pipeline end-to-end."""
    config = PipelineConfig(
        pam_window=pam_window,
        top_n_sgRNAs=top_n,
        homology_arm_length=arm_length,
        max_shift=max_shift,
        max_mismatches=max_mismatches,
        primer_tm_min=tm_min,
        primer_tm_max=tm_max,
        avoid_enzymes=avoid_enzymes or [],
    )

    if not force_fallback:
        config = check_external_tools(config)
    else:
        # --no-blast forces the k-mer fallback for sgRNA specificity only.
        # It must NOT disable RNAfold, which is used for secondary structure
        # prediction in primer/arm design and is independent of BLAST.
        config = check_external_tools(config)
        config.blast_available = False

    _log(f"Starting TAR-CRISPR pipeline...", verbose)
    _log(f"BLAST available: {config.blast_available}", verbose)
    _log(f"RNAfold available: {config.rnafold_available}", verbose)

    # Step 1: Ingestion and validation
    _log("Step 1: Reading input sequences...", verbose)
    bgc_record = read_sequence(genbank if genbank else bgc)
    vector_record = read_sequence(vector)

    genome_record = None
    genome_offsets = None
    if genome:
        genome_record, genome_offsets = read_genome(genome)
        if len(genome_offsets) > 1:
            _log(f"Genome: {len(genome_offsets)} scaffolds, "
                f"{len(genome_record.seq)} bp flattened (60N separators)", verbose)

    bgc_stats = validate_sequence(bgc_record)
    vector_stats = validate_sequence(vector_record)
    genome_stats = validate_sequence(genome_record) if genome_record else bgc_stats

    if not bgc_stats["valid"]:
        sys.exit(f"ERROR: BGC sequence has invalid characters: {bgc_stats['invalid_chars']}")
    if not vector_stats["valid"]:
        sys.exit(f"ERROR: Vector sequence has invalid characters: {vector_stats['invalid_chars']}")

    _log(f"BGC: {bgc_stats['length']} bp, GC={bgc_stats['gc_percent']}%", verbose)
    _log(f"Vector: {vector_stats['length']} bp, GC={vector_stats['gc_percent']}%", verbose)

    # Determine cluster bounds.
    # If the user gave a FASTA --bgc (no features) and omitted --start/--end,
    # extract_cluster_bounds cannot infer anything. Fail early with a clear
    # message instead of raising a ValueError deep in the call stack.
    has_features = bool(getattr(bgc_record, "features", None))

    # If the user has a genome AND a BGC record but did not say where in
    # the genome it is, locate the BGC via BLAST. This applies whether or
    # not the BGC has features:
    #   - no features -> the whole BGC is the cluster, its offset in the
    #     flat genome is the cluster location.
    #   - has features (antiSMASH) -> the features define the cluster
    #     RELATIVE TO THE BGC; we add the BGC's flat offset to translate.
    if start is None and end is None and genome_record is not None:
        _log("No --start/--end given; locating BGC in genome via BLAST...", verbose)
        try:
            scaffold_id, local_start, local_end, identity, locate_warnings = \
                locate_bgc_in_genome(bgc_record, genome_record, genome_offsets)
        except ValueError as e:
            sys.exit(
                f"ERROR: Could not locate the BGC in the given genome: {e}\n"
                "  - Pass --start/--end manually, or use --genbank with an "
                "annotated file."
            )
        for w in locate_warnings:
            _log(f"WARNING: {w}", verbose)
        _log(f"BGC located: {scaffold_id}, identity={identity:.1f}%", verbose)

        # Offset of the BGC in the flat genome
        bgc_flat_start, _ = translate_to_flat(scaffold_id, local_start, local_end,
                                              genome_offsets)

        # Determine the cluster bounds in the BGC's own coordinate space.
        # If the BGC has features (antiSMASH region, etc.), use them; if it
        # is a bare FASTA, the whole BGC is the cluster.
        if has_features:
            try:
                cluster_local = extract_cluster_bounds(bgc_record,
                                                       gene_kinds=gene_kinds)
            except ValueError as e:
                sys.exit(f"ERROR: Could not determine cluster bounds: {e}")
            local_start_c, local_end_c = cluster_local.start, cluster_local.end
        else:
            local_start_c, local_end_c = 0, len(bgc_record.seq)

        # Translate to flat genome coordinates
        start = bgc_flat_start + local_start_c
        end = bgc_flat_start + local_end_c
        _log(f"Cluster (flat): {start}..{end} "
             f"(local {local_start_c}..{local_end_c})", verbose)
        cluster_target = genome_record
    else:
        # Unchanged: start/end (explicit or from --genbank features) are
        # always relative to bgc_record, never to genome_record.
        cluster_target = bgc_record

    if (start is None or end is None) and not has_features and genbank is None:
        sys.exit(
            "ERROR: Cannot determine cluster bounds.\n"
            "  - You provided --bgc as a FASTA without --start/--end, and the file\n"
            "    has no features to parse coordinates from.\n"
            "  - Either pass --start and --end (0-based, half-open), or use\n"
            "    --genbank with an annotated file, or provide --genome so the\n"
            "    BGC can be located automatically via BLAST."
        )

    try:
        cluster = extract_cluster_bounds(cluster_target, start, end, gene_kinds=gene_kinds)
    except ValueError as e:
        sys.exit(f"ERROR: Could not determine cluster bounds: {e}")
    _log(f"Cluster: {cluster.name}, start={cluster.start}, end={cluster.end}", verbose)

    # Get flanking sequences
    if genome_record:
        upstream, cluster_seq, downstream = get_flanking_sequence(
            genome_record, cluster, config.pam_window
        )
    else:
        upstream, cluster_seq, downstream = get_flanking_sequence(
            bgc_record, cluster, config.pam_window
        )

    if not upstream:
        _log("WARNING: Upstream flank is empty (cluster near start of sequence)", verbose)
    if not downstream:
        _log("WARNING: Downstream flank is empty (cluster near end of sequence)", verbose)

    # Step 2: sgRNA design
    _log("Step 2: Designing sgRNAs...", verbose)
    genome_seq = str(genome_record.seq) if genome_record else str(bgc_record.seq)
    sgRNAs = design_sgRNAs(
        cluster_seq, upstream, downstream, genome_seq, config
    )
    _log(f"Found {len(sgRNAs['left'])} left, {len(sgRNAs['right'])} right candidates", verbose)

    # Fail early if either side has no sgRNA. Without cuts at both
    # boundaries, the fragment cannot be excised and downstream steps
    # would hit an AttributeError on NoneType.
    missing = [side for side in ("left", "right") if not sgRNAs[side]]
    if missing:
        sides = " and ".join(missing)
        hint = ""
        if any(s in missing for s in ("left",)) and cluster.start == 0:
            hint = (
                "\n  Hint: the cluster touches the start of the genome "
                "(start=0), so there is no upstream flank to search for "
                "PAM sites. Trim the cluster or provide a genome with "
                "extra upstream sequence."
            )
        if any(s in missing for s in ("right",)) and cluster.end >= len(str(genome_record.seq)):
            hint += (
                "\n  Hint: the cluster touches the end of the genome "
                "(end=genome_length), so there is no downstream flank. "
                "Trim the cluster or extend the genome."
            )
        sys.exit(
            f"ERROR: No sgRNA candidates found at the {sides} boundary. "
            f"Cannot excise a fragment.{hint}"
        )

    # Select sgRNAs
    if auto_select and sgRNAs["left"] and sgRNAs["right"]:
        selected = {"left": sgRNAs["left"][0], "right": sgRNAs["right"][0]}
        _log(f"Auto-selected left: {selected['left'].protospacer}", verbose)
        _log(f"Auto-selected right: {selected['right'].protospacer}", verbose)
    else:
        selected = _select_sgRNAs_interactive(sgRNAs, verbose)

    # Step 3: Extract fragment
    _log("Step 3: Extracting Cas9 fragment...", verbose)
    full_seq = genome_seq
    try:
        fragment = extract_fragment(full_seq, selected["left"], selected["right"])
        _log(f"Fragment: {fragment.fragment_length} bp, GC={fragment.gc_percent}%", verbose)
    except ValueError as e:
        sys.exit(f"ERROR: Could not extract fragment: {e}")

    # Step 4: Homology arms
    _log("Step 4: Designing homology arms...", verbose)
    homology_arms = design_homology_arms(
        fragment.sequence, genome_seq, str(vector_record.seq), config
    )
    _log(f"Left arm: {homology_arms['left'].sequence[:20]}...", verbose)
    _log(f"Right arm: {homology_arms['right'].sequence[-20:]}...", verbose)

    # Step 5: Tailed primers
    _log("Step 5: Designing tailed primers...", verbose)
    vec_seq_str = str(vector_record.seq).upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if vector_cut_left is None or vector_cut_right is None:
        cut_left, cut_right = _auto_detect_cut_site(vec_seq_str, vector_enzyme)
        if cut_left is not None and cut_right is not None:
            if vector_cut_left is None:
                vector_cut_left = cut_left
            if vector_cut_right is None:
                vector_cut_right = cut_right
    if vector_cut_right is None:
        vector_cut_right = vector_stats["length"]
    if vector_cut_left is None:
        vector_cut_left = 0
    _log(f"Vector cut: left={vector_cut_left}, right={vector_cut_right}", verbose)
    primers = design_tailed_primers(
        homology_arms["left"], homology_arms["right"],
        str(vector_record.seq), vector_cut_left, vector_cut_right, config
    )
    _log(f"Left primer Tm: {primers['left'].tm}C", verbose)
    _log(f"Right primer Tm: {primers['right'].tm}C", verbose)

    # Step 6: Assembly simulation
    _log("Step 6: Simulating assembly...", verbose)
    assembly_result = simulate_pydna_assembly(
        fragment.sequence, str(vector_record.seq),
        homology_arms["left"], homology_arms["right"],
        vector_cut_left, vector_cut_right,
    )
    _log(f"Assembly {'passed' if assembly_result.success else 'FAILED'}: "
         f"{assembly_result.final_size} bp, circular={assembly_result.circular}", verbose)

    # Step 7: Generate report
    _log("Step 7: Generating report...", verbose)
    Path(output).mkdir(parents=True, exist_ok=True)
    report_path = generate_report(
        sgRNAs=sgRNAs,
        fragment_seq=fragment.sequence,
        homology_arms=homology_arms,
        primers=primers,
        assembly_result=assembly_result,
        cluster_info=cluster,
        config=config,
        genome_stats=bgc_stats,
        output_dir=output,
        selected_sgRNAs=selected,
    )

    # Optional: render a circular map of the final construct.
    construct_png = Path(output) / "construct_map.png"
    try:
        generate_construct_map(
            final_sequence=assembly_result.final_sequence,
            vector_record=vector_record,
            bgc_record=bgc_record,
            vector_cut_left=vector_cut_left,
            vector_cut_right=vector_cut_right,
            bgc_len=len(fragment.sequence),
            output_png=str(construct_png),
            bgc_name=cluster.name or "BGC",
        )
    except Exception as e:
        _log(f"WARNING: Could not render construct map: {e}", verbose)
        construct_png = None

    _log(f"\nPipeline complete!", verbose)
    _log(f"Report: {report_path}", verbose)
    _log(f"SVG map: {Path(output) / 'cluster_map.svg'}", verbose)
    _log(f"Primers CSV: {Path(output) / 'primers.csv'}", verbose)
    if construct_png is not None:
        _log(f"Construct map: {construct_png}", verbose)

    if not assembly_result.success:
        _log("WARNING: Assembly simulation reported issues. Check report for details.", verbose)


def _log(msg: str, verbose: bool):
    """Print message if verbose is enabled."""
    if verbose:
        print(f"[TAR-CRISPR] {msg}", flush=True)


def _select_sgRNAs_interactive(sgRNAs: dict, verbose: bool) -> dict:
    """Interactive sgRNA selection (prompts user to pick from ranked list)."""
    selected = {}
    for end in ("left", "right"):
        end_name = "upstream (left)" if end == "left" else "downstream (right)"
        print(f"\n=== {end_name} sgRNA candidates (ranked) ===")
        candidates = sgRNAs[end]
        for i, sg in enumerate(candidates):
            print(f"  [{i+1}] Protospacer: {sg.protospacer}")
            print(f"      Strand: {sg.strand}, Cut: {sg.cut_position}, "
                  f"GC: {sg.gc_percent}%, Poly-T: {'Yes' if sg.polyt_flag else 'No'}")
        if candidates:
            choice = input(f"Select [{1}-{len(candidates)}] (default 1): ").strip()
            idx = int(choice) - 1 if choice else 0
            idx = max(0, min(idx, len(candidates) - 1))
            selected[end] = candidates[idx]
        else:
            print(f"  No candidates found for {end_name}!")
            selected[end] = None
    return selected


def _auto_detect_cut_site(vector_seq: str, enzyme_name: str) -> tuple:
    """Find a unique restriction site in the vector and return cut coordinates.

    Parameters
    ----------
    vector_seq : str
        The vector sequence (uppercased).
    enzyme_name : str
        Name of the restriction enzyme (e.g. "EcoRI").

    Returns
    -------
    tuple of (cut_left, cut_right) or (None, None) if not found.
    """
    try:
        from Bio.Restriction import Restriction
        from Bio.Seq import Seq
        enz = getattr(Restriction, enzyme_name)
        # Get the recognition site and cut position
        site = str(enz.site)
        # EcoRI cuts between G and AATTC (5'->3'), so cut at the boundary
        # The restriction enzyme site length
        site_len = len(site)
        # Find all occurrences
        positions = []
        for i in range(len(vector_seq) - site_len + 1):
            if vector_seq[i:i + site_len] == site:
                positions.append(i)
        if len(positions) == 0:
            return None, None
        if len(positions) == 1:
            pos = positions[0]
            # Cut at the enzyme's cut site within the recognition sequence
            # For EcoRI: G^AATTC, cut at position 1 within the site
            cut_pos = pos  # Cut at the start of the site (simplified: cut the whole site out)
            cut_left = pos
            cut_right = pos + site_len
            return cut_left, cut_right
        else:
            # Use the first unique site
            pos = positions[0]
            return pos, pos + site_len
    except (AttributeError, ImportError):
        return None, None


@app.command("sgrnas")
def sgrna_command(
    genome: str = typer.Option(..., "--genome", "-g", help="Genome FASTA"),
    start: int = typer.Option(0, "--start", help="BGC start (0-based)"),
    end: int = typer.Option(1000, "--end", help="BGC end (0-based)"),
    pam_window: int = typer.Option(500, "--pam-window", help="Flanking window"),
    top_n: int = typer.Option(5, "--top-n", help="Top N candidates"),
    max_mismatches: int = typer.Option(3, "--max-mismatches", help="Max mismatches"),
):
    """Quick sgRNA design for a given genomic region."""
    record = read_sequence(genome)
    cluster = extract_cluster_bounds(record, start, end)
    upstream, cluster_seq, downstream = get_flanking_sequence(record, cluster, pam_window)
    config = PipelineConfig(
        pam_window=pam_window,
        top_n_sgRNAs=top_n,
        max_mismatches=max_mismatches,
    )
    config = check_external_tools(config)
    genome_seq = str(record.seq)
    sgRNAs = design_sgRNAs(cluster_seq, upstream, downstream, genome_seq, config)

    print(f"\n=== Left (upstream) sgRNAs ===")
    for i, sg in enumerate(sgRNAs["left"]):
        print(f"  {i+1}. {sg.protospacer} | {sg.strand} strand | "
              f"cut={sg.cut_position} | GC={sg.gc_percent}% | Poly-T={'Y' if sg.polyt_flag else 'N'}")

    print(f"\n=== Right (downstream) sgRNAs ===")
    for i, sg in enumerate(sgRNAs["right"]):
        print(f"  {i+1}. {sg.protospacer} | {sg.strand} strand | "
              f"cut={sg.cut_position} | GC={sg.gc_percent}% | Poly-T={'Y' if sg.polyt_flag else 'N'}")


if __name__ == "__main__":
    app()
