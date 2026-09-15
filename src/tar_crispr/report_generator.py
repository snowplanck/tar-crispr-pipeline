"""Step 7: Report generation (Markdown + SVG visualization + CSV export)."""
import csv
import os
import io
from dataclasses import asdict
from typing import Optional, Any
from datetime import datetime

from tar_crispr.config import (
    PipelineConfig, ClusterInfo, PAMCandidate, HomologyArm,
    TailedPrimer, AssemblyResult,
)


def _format_sgRNA_table(sgRNAs: dict) -> str:
    """Format sgRNA candidates into a Markdown table."""
    lines = [
        "| End | Strand | Position | Protospacer (20nt) | PAM | Cut Pos | GC% | Poly-T |",
        "|------|--------|----------|-------------------|-----|---------|-----|--------|",
    ]
    for end_label, candidates in sgRNAs.items():
        end_name = {"left": "Upstream", "right": "Downstream"}.get(end_label, end_label)
        for sg in candidates:
            lines.append(
                f"| {end_name} | {sg.strand} | {sg.position} | "
                f"`{sg.protospacer}` | {sg.pam} | {sg.cut_position} | "
                f"{sg.gc_percent} | {'Yes' if sg.polyt_flag else 'No'} |"
            )
    return "\n".join(lines)


def _format_arm_table(arms: dict) -> str:
    """Format homology arms into a Markdown table."""
    lines = [
        "| End | Length | Sequence | GC% | Start | End | Unique | Structure | Issues |",
        "|------|--------|----------|-----|-------|-----|--------|-----------|--------|",
    ]
    for end_label, arm in arms.items():
        end_name = {"left": "Left (5')", "right": "Right (3')"}.get(end_label, end_label)
        issues_str = "; ".join(arm.issues) if arm.issues else "None"
        lines.append(
            f"| {end_name} | {arm.length} | `{arm.sequence}` | "
            f"{arm.gc_percent} | {arm.start_pos} | {arm.end_pos} | "
            f"{'Yes' if arm.uniqueness else 'No'} | "
            f"{'Yes' if arm.secondary_structure else 'No'} | {issues_str} |"
        )
    return "\n".join(lines)


def _format_primer_table(primers: dict) -> str:
    """Format tailed primers into a Markdown table."""
    lines = [
        "| Name | Sequence | Tail Len | Anneal Len | Tm (°C) | GC% | "
        "Hairpin | Self-Dimer | Cross-Dimer | Valid | Issues | Warnings |",
        "|------|----------|----------|------------|---------|-----|"
        "---------|------------|-------------|-------|--------|----------|",
    ]
    for key in ("left", "right"):
        p = primers[key]
        issues_str = "; ".join(p.issues) if p.issues else "—"
        warnings_str = "; ".join(getattr(p, "warnings", []) or []) or "—"
        lines.append(
            f"| {p.name} | `{p.sequence}` | {len(p.tail)} | "
            f"{len(p.annealing_region)} | "
            f"{p.tm} | {p.gc_percent} | {'Yes' if p.hairpin else 'No'} | "
            f"{'Yes' if p.self_dimer else 'No'} | "
            f"{'Yes' if p.cross_dimer else 'No'} | "
            f"{'Yes' if p.valid else 'No'} | {issues_str} | {warnings_str} |"
        )
    return "\n".join(lines)


def generate_svg_cluster_map(cluster_seq: str,
                             left_sg: PAMCandidate,
                             right_sg: PAMCandidate,
                             left_arm: HomologyArm,
                             right_arm: HomologyArm,
                             cluster_info: ClusterInfo) -> str:
    """Generate an SVG schematic of the BGC with cut sites and homology arms.

    Parameters
    ----------
    cluster_seq : str
        The BGC sequence.
    left_sg : PAMCandidate
        Left sgRNA.
    right_sg : PAMCandidate
        Right sgRNA.
    left_arm : HomologyArm
        Left homology arm.
    right_arm : HomologyArm
        Right homology arm.
    cluster_info : ClusterInfo
        Cluster coordinates in the genome.

    Returns
    -------
    str
        SVG string.
    """
    frag_len = len(cluster_seq)
    # Scale to a fixed width
    width = 800
    height = 200
    scale = width / max(frag_len, 1)

    # Coordinates (relative to fragment)
    left_cut_x = 0  # Left arm starts at 0
    right_cut_x = frag_len  # Right arm ends at frag_len

    arm_left_end = left_arm.length  # End of left arm
    arm_right_start = frag_len - right_arm.length  # Start of right arm

    def scale_x(pos: int) -> float:
        return pos * scale

    elements = []

    # Background
    elements.append(f'<rect x="0" y="0" width="{width}" height="{height}" '
                     f'fill="#f5f5f5" stroke="#333" stroke-width="1"/>')

    # BGC backbone (the full fragment)
    y_center = height / 2
    bar_y = y_center - 15
    bar_height = 30
    elements.append(f'<rect x="0" y="{bar_y}" width="{width}" '
                     f'height="{bar_height}" fill="none" stroke="#999" '
                     f'stroke-width="1" stroke-dasharray="4,2"/>')

    # Left homology arm (highlighted region at start)
    arm_l_w = scale_x(left_arm.length)
    elements.append(f'<rect x="0" y="{bar_y}" width="{arm_l_w}" '
                     f'height="{bar_height}" fill="#4CAF50" fill-opacity="0.3" '
                     f'stroke="#4CAF50" stroke-width="1"/>')
    elements.append(f'<text x="5" y="{bar_y - 5}" '
                     f'font-size="12" text-anchor="start" fill="#333">'
                     f'Left Arm ({left_arm.length}bp)</text>')

    # Right homology arm
    arm_r_x = scale_x(arm_right_start)
    arm_r_w = width - arm_r_x
    elements.append(f'<rect x="{arm_r_x}" y="{bar_y}" width="{arm_r_w}" '
                     f'height="{bar_height}" fill="#2196F3" fill-opacity="0.3" '
                     f'stroke="#2196F3" stroke-width="1"/>')
    elements.append(f'<text x="{width - 5}" '
                     f'y="{bar_y - 5}" font-size="12" text-anchor="end" '
                     f'fill="#333">Right Arm ({right_arm.length}bp)</text>')

    # Cas9 cut site markers
    cut_l_x = scale_x(left_arm.length)
    elements.append(f'<line x1="{cut_l_x}" y1="{bar_y}" x2="{cut_l_x}" '
                     f'y2="{bar_y + bar_height}" stroke="#F44336" stroke-width="2"/>')
    elements.append(f'<text x="{cut_l_x + 5}" y="{y_center}" font-size="12" '
                     f'fill="#F44336">Cas9 Left Cut</text>')

    cut_r_x = scale_x(arm_right_start)
    elements.append(f'<line x1="{cut_r_x}" y1="{bar_y}" x2="{cut_r_x}" '
                     f'y2="{bar_y + bar_height}" stroke="#F44336" stroke-width="2"/>')
    elements.append(f'<text x="{cut_r_x - 5}" y="{y_center}" font-size="12" '
                     f'text-anchor="end" fill="#F44336">Cas9 Right Cut</text>')

    # BGC core label
    bgc_x = scale_x(frag_len / 2)
    elements.append(f'<text x="{bgc_x}" y="{y_center + 30}" font-size="14" '
                     f'text-anchor="middle" fill="#333" font-weight="bold">'
                     f'BGC Core ({frag_len} bp)</text>')

    # Scale bar
    bar_x = 20
    bar_y2 = height - 30
    elements.append(f'<line x1="{bar_x}" y1="{bar_y2}" x2="{bar_x + 100}" '
                    f'y2="{bar_y2}" stroke="#333" stroke-width="2"/>')
    elements.append(f'<text x="{bar_x + 50}" y="{bar_y2 - 5}" font-size="11" '
                    f'text-anchor="middle" fill="#333">100 bp</text>')

    # Legend
    legend_y = height - 10
    elements.append(f'<rect x="20" y="{legend_y - 10}" width="12" height="12" '
                     f'fill="#4CAF50" fill-opacity="0.3" stroke="#4CAF50"/>')
    elements.append(f'<text x="38" y="{legend_y}" font-size="11" fill="#333">'
                     f'Left Homology Arm</text>')
    elements.append(f'<rect x="160" y="{legend_y - 10}" width="12" height="12" '
                     f'fill="#2196F3" fill-opacity="0.3" stroke="#2196F3"/>')
    elements.append(f'<text x="178" y="{legend_y}" font-size="11" fill="#333">'
                     f'Right Homology Arm</text>')
    elements.append(f'<line x1="340" y1="{legend_y - 4}" x2="352" y2="{legend_y - 4}" '
                    f'stroke="#F44336" stroke-width="2"/>')
    elements.append(f'<text x="358" y="{legend_y}" font-size="11" fill="#333">'
                     f'Cas9 Cut Site</text>')

    svg = (
        f'<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        + "\n".join(elements) + "\n</svg>"
    )
    return svg


def generate_report(sgRNAs: dict,
                    fragment_seq: str,
                    homology_arms: dict,
                    primers: dict,
                    assembly_result: AssemblyResult,
                    cluster_info: ClusterInfo,
                    config: PipelineConfig,
                    genome_stats: dict,
                    output_dir: str,
                    selected_sgRNAs: Optional[dict] = None) -> str:
    """Generate the final report in Markdown format with embedded SVG.

    Parameters
    ----------
    sgRNAs : dict
        Ranked sgRNA candidates for both ends (from design_sgRNAs).
    fragment_seq : str
        The Cas9-cut fragment sequence.
    homology_arms : dict
        Designed homology arms (from design_homology_arms).
    primers : dict
        Designed tailed primers (from design_tailed_primers).
    assembly_result : AssemblyResult
        Result from the assembly simulation.
    cluster_info : ClusterInfo
        Cluster coordinates.
    config : PipelineConfig
        Pipeline configuration used.
    genome_stats : dict
        Sequence validation stats from Step 1.
    output_dir : str
        Directory to write report files.
    selected_sgRNAs : dict, optional
        User-selected sgRNAs (overrides auto-selection).

    Returns
    -------
    str
        Path to the generated Markdown report.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Select sgRNAs (use user-selected or top-ranked)
    left_sg = selected_sgRNAs["left"] if selected_sgRNAs and "left" in selected_sgRNAs else sgRNAs["left"][0]
    right_sg = selected_sgRNAs["right"] if selected_sgRNAs and "right" in selected_sgRNAs else sgRNAs["right"][0]

    # Generate SVG cluster map
    svg_content = generate_svg_cluster_map(
        fragment_seq, left_sg, right_sg,
        homology_arms["left"], homology_arms["right"],
        cluster_info
    )
    svg_path = os.path.join(output_dir, "cluster_map.svg")
    with open(svg_path, "w", encoding="utf-8") as fh:
        fh.write(svg_content)

    # Export primers CSV
    from tar_crispr.primer_design import export_primers
    csv_path = os.path.join(output_dir, "primers.csv")
    export_primers(primers, csv_path)

    # Build Markdown report
    md = []
    md.append(f"# TAR-CRISPR Pipeline Report")
    md.append(f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md.append(f"\n**Pipeline version:** `{__import__('tar_crispr').__version__}`")
    md.append("\n---\n")

    # Section 1: Input Summary
    md.append("## 1. Input Sequence Summary\n")
    md.append(f"| Metric | Value |")
    md.append(f"|--------|-------|")
    md.append(f"| Sequence ID | {genome_stats.get('seqrecord_id', 'N/A')} |")
    md.append(f"| Description | {genome_stats.get('description', 'N/A')} |")
    md.append(f"| Total length | {genome_stats.get('length', 'N/A')} bp |")
    md.append(f"| GC content | {genome_stats.get('gc_percent', 'N/A')}% |")
    md.append(f"| Valid characters | {'Yes' if genome_stats.get('valid', False) else 'No'} |")
    if genome_stats.get("invalid_chars"):
        md.append(f"| Invalid characters | {', '.join(genome_stats['invalid_chars'])} |")
    else:
        md.append(f"| Invalid characters | None |")
    md.append(f"\n### Cluster Coordinates")
    md.append(f"- **Cluster:** {cluster_info.name}")
    md.append(f"- **Start (0-based):** {cluster_info.start}")
    md.append(f"- **End (0-based):** {cluster_info.end}")
    md.append(f"- **Size:** {cluster_info.end - cluster_info.start} bp")
    md.append(f"- **Strand:** {cluster_info.strand}")
    md.append("")

    # Section 2: sgRNA Design
    md.append("## 2. sgRNA Design (SpCas9 NGG)\n")
    md.append(f"**PAM window:** ±{config.pam_window} bp | "
              f"**Protospacer length:** {config.protospacer_length} nt | "
              f"**Top N:** {config.top_n_sgRNAs}\n")
    md.append(f"**Selected left sgRNA:** `{left_sg.protospacer}` ({left_sg.strand} strand, "
              f"cut at {left_sg.cut_position})\n")
    md.append(f"**Selected right sgRNA:** `{right_sg.protospacer}` ({right_sg.strand} strand, "
              f"cut at {right_sg.cut_position})\n")
    md.append(f"\n### All Candidate sgRNAs\n")
    md.append(_format_sgRNA_table(sgRNAs))
    md.append("")

    # Section 3: Fragment Definition
    md.append("## 3. Excised Fragment\n")
    md.append(f"| Property | Value |")
    md.append(f"|----------|-------|")
    md.append(f"| Fragment length | {len(fragment_seq)} bp |")
    md.append(f"| Fragment GC% | {round(_gc_simple(fragment_seq) * 100, 2)}% |")
    md.append(f"| Left cut position | {left_sg.cut_position} |")
    md.append(f"| Right cut position | {right_sg.cut_position} |")
    md.append(f"\n**Fragment sequence (first 100bp):**")
    md.append(f"```\n{fragment_seq[:100]}...\n```")
    md.append("")

    # Section 4: Homology Arms
    md.append("## 4. Homology Arm Design\n")
    md.append(f"**Arm length:** {config.homology_arm_length} bp | "
              f"**Max shift:** {config.max_shift} bp | "
              f"**Shift increment:** {config.shift_increment} bp\n")
    md.append(_format_arm_table(homology_arms))
    md.append("")

    # Section 5: Tailed Primers
    md.append("## 5. Tailed Primer Design\n")
    md.append(f"**Annealing region:** {config.primer_anneal_min}-{config.primer_anneal_max} nt | "
              f"**Tm range:** {config.primer_tm_min}-{config.primer_tm_max}°C\n")
    md.append(_format_primer_table(primers))
    md.append(f"\n**Primer CSV exported to:** `primers.csv`")
    md.append("")

    # Section 6: Assembly Simulation
    md.append("## 6. In Silico Assembly Simulation\n")
    md.append(f"| Property | Value |")
    md.append(f"|----------|-------|")
    md.append(f"| Assembly success | {'**PASSED**' if assembly_result.success else '**FAILED**'} |")
    md.append(f"| Final construct size | {assembly_result.final_size} bp |")
    md.append(f"| Circular | {'Yes' if assembly_result.circular else 'No'} |")
    md.append(f"| Final sequence length | {len(assembly_result.final_sequence)} bp |")
    if assembly_result.issues:
        md.append(f"\n### Assembly Issues")
        for issue in assembly_result.issues:
            md.append(f"- {issue}")
    else:
        md.append(f"\nNo issues detected. Assembly verified.")
    md.append("")

    # Section 7: Visualization
    md.append("## 7. BGC Schematic Map\n")
    md.append(f"![Cluster Map](cluster_map.svg)\n")
    md.append(f"*Green: left homology arm | Blue: right homology arm | "
              f"Red lines: Cas9 cut sites*\n")

    # Section 8: Configuration
    md.append("\n## 8. Pipeline Configuration\n")
    md.append(f"```")
    for key, value in asdict(config).items():
        if value is not None:
            md.append(f"{key}: {value}")
    md.append(f"```")

    # Footer
    md.append("\n---\n")
    md.append(f"*Generated by tar-crispr-pipeline v{__import__('tar_crispr').__version__}*\n")

    # Write report
    report_path = os.path.join(output_dir, "report.md")
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))

    # Also write HTML version
    html_path = os.path.join(output_dir, "report.html")
    _write_html_report(report_path, html_path, svg_path, csv_path, md)

    return report_path


def _gc_simple(seq: str) -> float:
    """Simple GC fraction."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return gc / len(seq)


def _write_html_report(md_path: str, html_path: str,
                       svg_path: str, csv_path: str, md_lines: list[str]) -> None:
    """Write an HTML version of the report with embedded SVG."""
    # Read the SVG
    with open(svg_path, "r", encoding="utf-8") as fh:
        svg_content = fh.read()

    # Convert markdown tables to HTML (simple approach)
    html_body = "\n".join(md_lines)

    # Replace markdown headers
    for i in range(1, 7):
        prefix = "#" * i
        html_body = html_body.replace(
            f"\n{prefix} ",
            f"\n<h{i}>"
        )
        # Close tag - tricky with simple replacement, use a different approach

    # Just wrap in basic HTML with the SVG inline
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TAR-CRISPR Pipeline Report</title>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; }}
  h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 10px; }}
  h2 {{ color: #2471a3; margin-top: 30px; }}
  table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
  th {{ background-color: #d6eaf8; }}
  tr:nth-child(even) {{ background-color: #f9f9f9; }}
  code {{ background: #f4f4f4; padding: 2px 6px; border-radius: 3px; }}
  .passed {{ color: green; font-weight: bold; }}
  .failed {{ color: red; font-weight: bold; }}
  .svg-container {{ text-align: center; margin: 20px 0; }}
  .footer {{ margin-top: 40px; padding-top: 10px; border-top: 1px solid #ddd; color: #888; font-size: 12px; }}
</style>
</head>
<body>
<h1>TAR-CRISPR Pipeline Report</h1>
<p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
{svg_content}
<hr>
<div class="footer">
<p>Generated by tar-crispr-pipeline. Primer CSV: <code>{csv_path}</code></p>
</div>
</body>
</html>"""

    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)


__all__ = [
    "generate_report",
    "generate_svg_cluster_map",
]
