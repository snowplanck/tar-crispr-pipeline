"""Tests for Cas9 cut-site scanning, modes, and the coordinate regressions."""
import random

import pytest

from tar_crispr.config import PipelineConfig, PAMCandidate
from tar_crispr.cut_specificity import (
    build_site_index, find_cut_sites, validate_guide_pair, select_valid_pair,
    build_index_from_fasta,
)
from tar_crispr.fragment_ends import extract_fragment
from tar_crispr.pam_finder import (
    design_sgRNAs, find_pam_sites, reverse_complement,
)

PROTO = "ACGTTGCAAGCTTGCAGTCA"


def rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def cand(proto, cut, strand="+"):
    return PAMCandidate(position=cut, strand=strand, protospacer=proto, pam="AGG",
                        cut_position=cut, gc_percent=50.0, polyt_flag=False)


# ---------------------------------------------------------------- cut coords
class TestCutCoordinates:
    def test_forward_cut_is_3bp_upstream_of_pam(self):
        seq = "TTTTTT" + PROTO + "AGG" + "TTTTTT"
        fwd = [c for c in find_pam_sites(seq) if c.strand == "+" and c.protospacer == PROTO]
        assert fwd[0].cut_position == 6 + 17

    def test_reverse_cut_is_3bp_from_pam_inside_protospacer(self):
        # Regression: used to be reported at the FAR end of the protospacer (+17).
        top = "TTTTTT" + reverse_complement("AGG") + reverse_complement(PROTO) + "TTTTTT"
        rev = [c for c in find_pam_sites(top) if c.strand == "-" and c.protospacer == PROTO]
        assert rev[0].cut_position == 6 + 3 + 3

    def test_index_and_finder_agree_on_both_strands(self):
        fwd = "TTTTTT" + PROTO + "AGG" + "TTTTTT"
        rev = "TTTTTT" + reverse_complement("AGG") + reverse_complement(PROTO) + "TTTTTT"
        for seq, strand in ((fwd, "+"), (rev, "-")):
            idx = build_site_index(seq)
            hit = [h for h in find_cut_sites(PROTO, idx) if h.strand == strand]
            finder = [c for c in find_pam_sites(seq) if c.strand == strand and c.protospacer == PROTO]
            assert hit and finder
            assert hit[0].cut_position == finder[0].cut_position

    def test_offset_applied(self):
        seq = "TTTTTT" + PROTO + "AGG" + "TTTTTT"
        h = find_cut_sites(PROTO, build_site_index(seq, offset=1000))
        assert h[0].cut_position == 1000 + 23


# ---------------------------------------------------------------- site rules
class TestSiteClassification:
    def _idx(self, site, pam="AGG", offset=0):
        return build_site_index("TTTTTT" + site + pam + "TTTTTT", offset=offset)

    def test_perfect_ngg_is_high(self):
        assert find_cut_sites(PROTO, self._idx(PROTO))[0].severity == "high"

    def test_distal_mismatches_within_limit_are_high(self):
        site = "T" + PROTO[1:]                      # 1 mismatch, far from PAM
        s = find_cut_sites(PROTO, self._idx(site))[0]
        assert (s.severity, s.mismatches, s.seed_mismatches) == ("high", 1, 0)

    def test_too_many_distal_mismatches_ignored(self):
        site = "TGCA" + PROTO[4:]      # 4 mismatches (A/C/G/T all differ)
        assert find_cut_sites(PROTO, self._idx(site), max_mismatches=3) == []

    def test_seed_mismatch_is_medium(self):
        site = PROTO[:-1] + ("A" if PROTO[-1] != "A" else "C")
        s = find_cut_sites(PROTO, self._idx(site))[0]
        assert (s.severity, s.seed_mismatches) == ("medium", 1)

    def test_nag_pam_is_medium(self):
        s = find_cut_sites(PROTO, self._idx(PROTO, pam="AAG"))[0]
        assert (s.pam_class, s.severity) == ("NAG", "medium")

    def test_no_pam_no_site(self):
        assert find_cut_sites(PROTO, self._idx(PROTO, pam="ATT")) == []

    def test_circular_finds_site_spanning_origin(self):
        plasmid = (PROTO[10:] + "AGG") + rnd(300, 1) + PROTO[:10]
        assert find_cut_sites(PROTO, build_site_index(plasmid, circular=True)) != []
        assert find_cut_sites(PROTO, build_site_index(plasmid, circular=False)) == []


# ---------------------------------------------------------------- design
def _scenario(cluster_start=8000, cluster_end=12000, seed=1):
    genome = rnd(20000, seed)
    cs, ce = cluster_start, cluster_end
    return genome, genome[cs - 500:cs], genome[cs:ce], genome[ce:ce + 500], cs, ce


class TestCoordinateRegression:
    def test_cuts_are_in_genome_coordinates(self):
        # Regression: cuts were local (0-based on upstream+cluster+downstream),
        # so extract_fragment sliced the wrong part of the genome.
        genome, up, cl, dn, cs, ce = _scenario()
        cfg = PipelineConfig(blast_available=False)
        sg = design_sgRNAs(cl, up, dn, genome, cfg, upstream_start=cs - len(up))
        L, R = sg["left"][0], sg["right"][0]
        assert cs - 500 <= L.cut_position <= cs
        assert ce <= R.cut_position <= ce + 500
        frag = extract_fragment(genome, L, R)
        assert cl in frag.sequence


class TestInternalCutFilter:
    def _plant(self):
        genome, up, cl, dn, cs, ce = _scenario(seed=7)
        cfg = PipelineConfig(blast_available=False, top_n_sgRNAs=50)
        base = design_sgRNAs(cl, up, dn, genome, cfg, upstream_start=cs - len(up))
        victim = base["left"][0]
        planted = cl[:1500] + victim.protospacer + "TGG" + cl[1523:]
        return victim, up, planted, dn, cs, ce, cfg

    def test_candidate_cutting_inside_bgc_is_rejected(self):
        victim, up, planted, dn, cs, ce, cfg = self._plant()
        genome = up + planted + dn
        rej = {}
        sg = design_sgRNAs(planted, up, dn, genome, cfg,
                           upstream_start=0, rejected_out=rej)
        assert all(c.protospacer != victim.protospacer for c in sg["left"])
        assert any(c.protospacer == victim.protospacer for c, _ in rej["left"])

    def test_override_keeps_candidate_with_warning(self):
        victim, up, planted, dn, cs, ce, cfg = self._plant()
        cfg.exclude_internal_cuts = False
        genome = up + planted + dn
        sg = design_sgRNAs(planted, up, dn, genome, cfg, upstream_start=0)
        kept = [c for c in sg["left"] if c.protospacer == victim.protospacer]
        assert kept and kept[0].internal_cuts
        assert any("additional cut" in w for w in kept[0].warnings)


# ---------------------------------------------------------------- modes
class TestModes:
    def _pair(self):
        genome = rnd(6000, 3)
        left = cand(genome[1000:1020], 1017)
        right = cand(genome[4980:5000], 4997)
        return genome, left, right

    def test_invalid_mode_raises(self):
        genome, l, r = self._pair()
        with pytest.raises(ValueError):
            validate_guide_pair(l, r, genome, PipelineConfig(mode="bogus"))

    def test_in_vitro_skips_vector_and_yeast(self):
        genome, l, r = self._pair()
        vec = "TTTT" + l.protospacer + "AGG" + "TTTT"        # would be cut in vivo
        rep = validate_guide_pair(l, r, genome, PipelineConfig(mode="in-vitro"), vector_seq=vec)
        assert rep.ok and any("in-vitro" in n for n in rep.notes)

    def test_in_vivo_vector_site_is_hard_issue(self):
        genome, l, r = self._pair()
        vec = rnd(200, 5) + l.protospacer + "AGG" + rnd(200, 6)
        rep = validate_guide_pair(l, r, genome, PipelineConfig(mode="in-vivo"), vector_seq=vec)
        assert not rep.ok and any("capture vector" in m for m in rep.hard_issues)

    def test_in_vivo_without_yeast_genome_warns(self):
        genome, l, r = self._pair()
        rep = validate_guide_pair(l, r, genome, PipelineConfig(mode="in-vivo"),
                                  vector_seq=rnd(300, 9))
        assert any("yeast genome not provided" in w for w in rep.warnings)

    def test_in_vivo_yeast_site_is_hard_issue(self, tmp_path):
        genome, l, r = self._pair()
        fa = tmp_path / "yeast.fa"
        fa.write_text(">chrI\n" + rnd(500, 11) + l.protospacer + "TGG" + rnd(500, 12) + "\n")
        yidx = build_index_from_fasta(str(fa))
        rep = validate_guide_pair(l, r, genome, PipelineConfig(mode="in-vivo"),
                                  vector_seq=rnd(300, 9), yeast_index=yidx)
        assert not rep.ok and any("yeast genome" in m for m in rep.hard_issues)

    def test_polyt_penalised_only_in_vivo(self):
        from tar_crispr.pam_finder import rank_sgRNAs
        polyt = PAMCandidate(0, "+", "GTTTTACGTACGTACGTACG", "AGG", 17, 45.0, True)
        clean = PAMCandidate(50, "+", "GACGTACGTACGTACGTACG", "AGG", 67, 55.0, False)
        def first(mode):
            ranked = rank_sgRNAs([dataclass_copy(polyt), dataclass_copy(clean)], "G" * 500,
                                 PipelineConfig(blast_available=False, mode=mode))
            return ranked[0].protospacer
        assert first("in-vivo") == clean.protospacer      # -15 for the TTTT run
        assert first("in-vitro") == polyt.protospacer     # tie, input order kept


def dataclass_copy(c):
    import dataclasses
    return dataclasses.replace(c, warnings=[], internal_cuts=[])


# ---------------------------------------------------------------- pair pick
class TestSelectValidPair:
    def test_skips_pair_with_internal_cut_and_picks_next(self):
        genome = rnd(6000, 21)
        l1 = cand(genome[1000:1020], 1017)
        l2 = cand(genome[1100:1120], 1117)
        r1 = cand(genome[4980:5000], 4997)
        # plant l1's target inside the fragment (perfect site + NGG)
        genome = genome[:3000] + l1.protospacer + "TGG" + genome[3023:]
        cfg = PipelineConfig(mode="in-vitro")
        pair, rep, tried = select_valid_pair({"left": [l1, l2], "right": [r1]}, genome, cfg)
        assert pair["left"] is l2 and rep.ok and tried == 2

    def test_returns_none_when_all_fail(self):
        genome = rnd(6000, 22)
        l1 = cand(genome[1000:1020], 1017)
        r1 = cand(genome[4980:5000], 4997)
        genome = genome[:3000] + l1.protospacer + "TGG" + genome[3023:]
        pair, rep, tried = select_valid_pair({"left": [l1], "right": [r1]}, genome,
                                             PipelineConfig(mode="in-vitro"))
        assert pair is None and rep is not None and not rep.ok


# ---------------------------------------------------------------- CLI e2e
class TestCliModes:
    """End-to-end: cuts must flank the cluster, and both modes must run."""

    def _run(self, tmp_path, mode):
        import re
        import subprocess
        import sys
        from pathlib import Path
        root = Path(__file__).resolve().parent.parent
        td = root / "test_data"
        out = tmp_path / mode
        proc = subprocess.run(
            [sys.executable, "-m", "tar_crispr.cli", "run",
             "--bgc", str(td / "synthetic_bgc.fasta"),
             "--vector", str(td / "synthetic_vector.fasta"),
             "--genome", str(td / "synthetic_genome.fasta"),
             "--genbank", str(td / "synthetic_bgc.gbk"),
             "--output", str(out), "--no-blast", "--mode", mode],
            capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        report = (out / "report.md").read_text()
        cuts = [int(x) for x in re.findall(r"cut at (\d+)", report)]
        return report, cuts

    @pytest.mark.parametrize("mode", ["in-vitro", "in-vivo"])
    def test_cuts_flank_cluster_and_report_has_safety_section(self, tmp_path, mode):
        report, cuts = self._run(tmp_path, mode)
        left, right = cuts[0], cuts[1]
        # synthetic cluster is 1000-4000 in genome coordinates
        assert left <= 1000 and right >= 4000
        assert f"Cas9 cut-site safety (mode: `{mode}`)" in report

    def test_invalid_mode_rejected(self, tmp_path):
        import subprocess
        import sys
        proc = subprocess.run([sys.executable, "-m", "tar_crispr.cli", "run",
                               "--bgc", "x", "--vector", "y", "--mode", "bogus"],
                              capture_output=True, text=True)
        assert proc.returncode != 0 and "--mode must be" in (proc.stdout + proc.stderr)
