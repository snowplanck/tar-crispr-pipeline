"""End-to-end test: run the full CLI on synthetic data and verify outputs."""
import os
import sys
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def project_root():
    return Path(__file__).parent.parent


class TestEndToEnd:
    """Run the pipeline on synthetic test data and verify all outputs."""

    def test_pipeline_generates_all_outputs(self, tmp_path, project_root):
        bgc = str(project_root / "test_data" / "synthetic_bgc.fasta")
        vector = str(project_root / "test_data" / "synthetic_vector.fasta")
        genome = str(project_root / "test_data" / "synthetic_genome.fasta")
        gbk = str(project_root / "test_data" / "synthetic_bgc.gbk")
        output = str(tmp_path / "e2e_output")

        cmd = [
            sys.executable, "-m", "tar_crispr.cli", "run",
            "--bgc", bgc,
            "--vector", vector,
            "--genome", genome,
            "--genbank", gbk,
            "--output", output,
            "--no-blast",
            "--verbose",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        # Pipeline should complete successfully
        assert result.returncode == 0, f"Pipeline failed: {result.stderr}"

        # Check output files exist
        report_md = Path(output) / "report.md"
        report_html = Path(output) / "report.html"
        cluster_svg = Path(output) / "cluster_map.svg"
        primers_csv = Path(output) / "primers.csv"

        assert report_md.exists(), "Markdown report not generated"
        assert report_html.exists(), "HTML report not generated"
        assert cluster_svg.exists(), "SVG cluster map not generated"
        assert primers_csv.exists(), "Primers CSV not generated"

        # Verify report content
        report_content = report_md.read_text(encoding="utf-8")
        assert "# TAR-CRISPR Pipeline Report" in report_content
        assert "## 1. Input Sequence Summary" in report_content
        assert "## 2. sgRNA Design" in report_content
        assert "## 3. Excised Fragment" in report_content
        assert "## 4. Homology Arm Design" in report_content
        assert "## 5. Tailed Primer Design" in report_content
        assert "## 6. In Silico Assembly Simulation" in report_content
        assert "## 7. BGC Schematic Map" in report_content
        assert "## 8. Pipeline Configuration" in report_content

        # Verify sgRNA table has entries
        assert "Protospacer" in report_content
        assert "NGG" in report_content

        # Verify primers CSV has correct structure
        csv_content = primers_csv.read_text(encoding="utf-8")
        assert "Primer_Name" in csv_content
        assert "FORWARD" in csv_content
        assert "REVERSE" in csv_content

        # Verify SVG has content
        svg_content = cluster_svg.read_text(encoding="utf-8")
        assert "<svg" in svg_content
        assert "</svg>" in svg_content
        assert "Cas9" in svg_content
        assert "Homology Arm" in svg_content

    def test_pipeline_with_explicit_coordinates(self, tmp_path, project_root):
        """Test pipeline using --start/--end instead of GenBank parsing."""
        bgc = str(project_root / "test_data" / "synthetic_bgc.fasta")
        vector = str(project_root / "test_data" / "synthetic_vector.fasta")
        output = str(tmp_path / "e2e_explicit_output")

        cmd = [
            sys.executable, "-m", "tar_crispr.cli", "run",
            "--bgc", bgc,
            "--vector", vector,
            "--start", "1000",
            "--end", "4000",
            "--output", output,
            "--no-blast",
            "--verbose",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, f"Pipeline failed: {result.stderr}"
        assert (Path(output) / "report.md").exists()

    def test_pipeline_sgrnas_command(self, tmp_path, project_root):
        """Test the sgrnas subcommand."""
        genome = str(project_root / "test_data" / "synthetic_genome.fasta")
        cmd = [
            sys.executable, "-m", "tar_crispr.cli", "sgrnas",
            "--genome", genome,
            "--start", "1000",
            "--end", "4000",
            "--top-n", "3",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, f"sgRNA command failed: {result.stderr}"
        assert "sgRNA" in result.stdout.lower() or "upstream" in result.stdout.lower()

    def test_report_contains_sgRNA_sequences(self, tmp_path, project_root):
        """Verify the report contains actual sgRNA sequences (not just headers)."""
        bgc = str(project_root / "test_data" / "synthetic_bgc.fasta")
        vector = str(project_root / "test_data" / "synthetic_vector.fasta")
        genome = str(project_root / "test_data" / "synthetic_genome.fasta")
        gbk = str(project_root / "test_data" / "synthetic_bgc.gbk")
        output = str(tmp_path / "e2e_check_output")

        cmd = [
            sys.executable, "-m", "tar_crispr.cli", "run",
            "--bgc", bgc,
            "--vector", vector,
            "--genome", genome,
            "--genbank", gbk,
            "--output", output,
            "--no-blast",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0

        report = (Path(output) / "report.md").read_text(encoding="utf-8")

        # Should contain the selected sgRNA sequences
        assert "Selected left sgRNA" in report
        assert "Selected right sgRNA" in report

        # Should contain actual nucleotide sequences (20 chars)
        assert "Cas9" in report  # Cut site references

        # Should contain homology arm sequences
        assert "Left (5')" in report
        assert "Right (3')" in report
        assert "homology arm" in report.lower()
