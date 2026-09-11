"""Tests for report_generator.py (Step 7)."""
import os
import pytest

from tar_crispr.report_generator import (
    generate_report,
    generate_svg_cluster_map,
    _format_sgRNA_table,
    _format_arm_table,
    _format_primer_table,
    _gc_simple,
)
from tar_crispr.config import (
    PipelineConfig, ClusterInfo, PAMCandidate, HomologyArm,
    TailedPrimer, AssemblyResult,
)


def make_sg(p, strand="+", cut_pos=100, gc=50.0, polyt=False):
    return PAMCandidate(
        position=p, strand=strand, protospacer="ACGTACGTACGTACGTACGT",
        pam="NGG", cut_position=cut_pos, gc_percent=gc, polyt_flag=polyt,
    )


def make_arm(seq, end="left"):
    return HomologyArm(
        sequence=seq, end=end, length=len(seq),
        gc_percent=50.0, start_pos=0, end_pos=len(seq),
        uniqueness=True, secondary_structure=False, issues=[],
    )


class TestGcSimple:
    def test_all_gc(self):
        assert _gc_simple("GCGCGC") == 1.0

    def test_no_gc(self):
        assert _gc_simple("AATATA") == 0.0

    def test_mixed(self):
        assert _gc_simple("GCAT") == 0.5

    def test_empty(self):
        assert _gc_simple("") == 0.0


class TestFormatTables:
    def test_sgRNA_table(self):
        sgRNAs = {"left": [make_sg(100)], "right": [make_sg(500, "right")]}
        table = _format_sgRNA_table(sgRNAs)
        assert "Upstream" in table
        assert "Downstream" in table
        assert "ACGTACGTACGTACGTACGT" in table

    def test_arm_table(self):
        arms = {"left": make_arm("GCGCGCGCATATATAT"), "right": make_arm("TATATATACGCGCGCG", "right")}
        table = _format_arm_table(arms)
        assert "Left" in table
        assert "Right" in table

    def test_primer_table(self):
        left_p = TailedPrimer(
            name="FORWARD_left", sequence="GCGCGCGCATATATATATACGTACGT",
            tail="GCGCGCGCATATATATAT", annealing_region="GTACGT",
            tm=60.0, gc_percent=50.0, hairpin=False, self_dimer=False, cross_dimer=False,
            issues=[], valid=True,
        )
        right_p = TailedPrimer(
            name="REVERSE_right", sequence="TATATATACGCGCGCGCGCGTACGT",
            tail="TATATATACGCGCGCGCGCG", annealing_region="GTACGT",
            tm=60.0, gc_percent=50.0, hairpin=False, self_dimer=False, cross_dimer=False,
            issues=[], valid=True,
        )
        primers = {"left": left_p, "right": right_p}
        table = _format_primer_table(primers)
        assert "FORWARD_left" in table
        assert "REVERSE_right" in table


class TestSvgClusterMap:
    def test_svg_generation(self):
        cluster_seq = "A" * 500 + "B" * 200 + "C" * 500
        left_sg = make_sg(100, cut_pos=50)
        right_sg = make_sg(800, "right", cut_pos=1150)
        left_arm = make_arm("GCGCGCGCATATATAT", "left")
        right_arm = make_arm("TATATATACGCGCGCG", "right")
        cluster_info = ClusterInfo(start=500, end=1500)
        svg = generate_svg_cluster_map(cluster_seq, left_sg, right_sg,
                                       left_arm, right_arm, cluster_info)
        assert "<?xml" in svg
        assert "<svg" in svg
        assert "</svg>" in svg
        assert "Cas9" in svg or "Cas9" in svg.lower()

    def test_svg_contains_arms(self):
        cluster_seq = "A" * 500 + "B" * 200 + "C" * 500
        left_sg = make_sg(100, cut_pos=50)
        right_sg = make_sg(800, "right", cut_pos=1150)
        left_arm = make_arm("GCGCGCGCATATATATAT", "left")
        right_arm = make_arm("TATATATATACGCGCGCGC", "right")
        cluster_info = ClusterInfo(start=500, end=1500)
        svg = generate_svg_cluster_map(cluster_seq, left_sg, right_sg,
                                       left_arm, right_arm, cluster_info)
        # Check for key labels
        assert "Homology Arm" in svg


class TestGenerateReport:
    def test_full_report(self, tmp_path):
        sgRNAs = {"left": [make_sg(100), make_sg(200)], "right": [make_sg(500, "right")]}
        fragment = "GCGCGCGCATATATATAT" + "BGC" * 100 + "TATATATACGCGCGCGC"
        arms = {"left": make_arm("GCGCGCGCATATATATAT"), "right": make_arm("TATATATACGCGCGCGC", "right")}
        left_p = TailedPrimer(
            name="FORWARD_left", sequence="GCGCGCGCATATATATATACGT",
            tail="GCGCGCGCATATATATAT", annealing_region="ACGT",
            tm=60.0, gc_percent=50.0, hairpin=False, self_dimer=False, cross_dimer=False,
            issues=[], valid=True,
        )
        right_p = TailedPrimer(
            name="REVERSE_right", sequence="TATATATACGCGCGCGCGCGT",
            tail="TATATATACGCGCGCGCGCG", annealing_region="GT",
            tm=60.0, gc_percent=50.0, hairpin=False, self_dimer=False, cross_dimer=False,
            issues=[], valid=True,
        )
        primers = {"left": left_p, "right": right_p}
        assembly = AssemblyResult(
            success=True, final_size=1000, circular=True,
            issues=[], final_sequence=fragment + "VECTOR_BACKBONE",
        )
        cluster_info = ClusterInfo(start=500, end=1500, name="testBGC")
        config = PipelineConfig(blast_available=False, rnafold_available=False)
        genome_stats = {"length": 6000, "gc_percent": 50.0, "valid": True,
                        "invalid_chars": [], "seqrecord_id": "test_genome"}

        report_path = generate_report(
            sgRNAs, fragment, arms, primers, assembly,
            cluster_info, config, genome_stats,
            str(tmp_path)
        )

        assert os.path.exists(report_path)
        with open(report_path) as fh:
            content = fh.read()

        assert "# TAR-CRISPR Pipeline Report" in content
        assert "## 1. Input Sequence Summary" in content
        assert "## 2. sgRNA Design" in content
        assert "## 3. Excised Fragment" in content
        assert "## 4. Homology Arm Design" in content
        assert "## 5. Tailed Primer Design" in content
        assert "## 6. In Silico Assembly Simulation" in content
        assert "## 7. BGC Schematic Map" in content
        assert "## 8. Pipeline Configuration" in content

        # Check SVG file exists
        assert os.path.exists(str(tmp_path / "cluster_map.svg"))

        # Check CSV file exists
        assert os.path.exists(str(tmp_path / "primers.csv"))

    def test_report_with_failed_assembly(self, tmp_path):
        sgRNAs = {"left": [make_sg(100)], "right": [make_sg(500, "right")]}
        fragment = "ACGT" * 100
        arms = {"left": make_arm("ACGTACGTACGTACGT"), "right": make_arm("ACGTACGTACGTACGT", "right")}
        left_p = TailedPrimer(
            name="FORWARD_left", sequence="ACGTACGTACGTACGTACGT",
            tail="ACGTACGTACGTACGT", annealing_region="ACGT",
            tm=60.0, gc_percent=50.0, hairpin=False, self_dimer=False, cross_dimer=False,
            issues=[], valid=True,
        )
        right_p = TailedPrimer(
            name="REVERSE_right", sequence="ACGTACGTACGTACGTACGT",
            tail="ACGTACGTACGTACGT", annealing_region="ACGT",
            tm=60.0, gc_percent=50.0, hairpin=False, self_dimer=False, cross_dimer=False,
            issues=[], valid=True,
        )
        primers = {"left": left_p, "right": right_p}
        assembly = AssemblyResult(
            success=False, final_size=0, circular=False,
            issues=["Assembly failed: junction mismatch"],
            final_sequence="",
        )
        cluster_info = ClusterInfo(start=100, end=500)
        config = PipelineConfig()
        genome_stats = {"length": 1000, "gc_percent": 50.0, "valid": True,
                        "invalid_chars": [], "seqrecord_id": "test"}

        report_path = generate_report(
            sgRNAs, fragment, arms, primers, assembly,
            cluster_info, config, genome_stats,
            str(tmp_path)
        )
        assert os.path.exists(report_path)
        with open(report_path) as fh:
            content = fh.read()
        assert "FAILED" in content
