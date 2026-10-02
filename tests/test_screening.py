"""Tests for colony-PCR screening primers and the relative arm-GC limits."""
import csv
import random

import primer3
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from tar_crispr.config import PipelineConfig
from tar_crispr.homology_arms import gc_limits, validate_gc
from tar_crispr.screening import (
    SIZE_LADDER, design_screening, export_screening, find_marker_gene,
    find_primer_sites, predict_amplicons,
)

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def rnd(n, seed, gc=False):
    r = random.Random(seed)
    w = [14, 36, 36, 14] if gc else [25, 25, 25, 25]
    return "".join(r.choices("ACGT", weights=w, k=n))


@pytest.fixture(scope="module")
def construct():
    v5, frag, v3 = rnd(3000, 1), rnd(37600, 2, gc=True), rnd(2500, 3)
    return {"v5": v5, "frag": frag, "v3": v3, "seq": v5 + frag + v3, "JL": len(v5)}


# ------------------------------------------------------------ in-silico PCR
class TestInSilicoPcr:
    SEQ = rnd(2000, 5)

    def test_forward_and_reverse_sites(self):
        f, r = self.SEQ[100:120], rc(self.SEQ[500:520])
        assert ("+", 100, 120, 0) in find_primer_sites(f, self.SEQ)
        assert ("-", 500, 520, 0) in find_primer_sites(r, self.SEQ)

    def test_product_size(self):
        f, r = self.SEQ[100:120], rc(self.SEQ[500:520])
        prods = predict_amplicons({"F": f, "R": r}, self.SEQ)
        assert [(p["forward"], p["reverse"], p["size"]) for p in prods] == [("F", "R", 420)]

    def test_mismatch_tolerance_and_exact_3prime_anchor(self):
        f = self.SEQ[100:120]
        near = f[:3] + ("A" if f[3] != "A" else "C") + f[4:]      # 1 mismatch, 5' side
        assert find_primer_sites(near, self.SEQ)[0][3] == 1
        bad3 = f[:-2] + ("A" if f[-2] != "A" else "C") + f[-1]    # mismatch inside the 3' 10-mer
        assert not [s for s in find_primer_sites(bad3, self.SEQ) if s[0] == "+" and s[1] == 100]

    def test_too_many_mismatches_ignored(self):
        f = self.SEQ[100:120]
        mut = "".join(("A" if b != "A" else "C") if i in (0, 2, 4, 6) else b for i, b in enumerate(f))
        assert find_primer_sites(mut, self.SEQ) == []

    def test_circular_product_spanning_origin(self):
        f, r = self.SEQ[1900:1920], rc(self.SEQ[100:120])
        prods = predict_amplicons({"F": f, "R": r}, self.SEQ, circular=True)
        assert prods and prods[0]["size"] == 220
        assert predict_amplicons({"F": f, "R": r}, self.SEQ, circular=False) == []

    def test_site_spanning_origin_found_only_when_circular(self):
        p = self.SEQ[-8:] + self.SEQ[:12]
        assert find_primer_sites(p, self.SEQ, circular=True)
        assert not find_primer_sites(p, self.SEQ, circular=False)

    def test_wrong_orientation_gives_no_product(self):
        f, r = self.SEQ[500:520], rc(self.SEQ[100:120])           # reverse primer upstream of forward
        prods = predict_amplicons({"F": f, "R": r}, self.SEQ, circular=False)
        assert prods == []


# ------------------------------------------------------------ panel design
class TestDesign:
    def _design(self, c, **kw):
        return design_screening(c["seq"], c["JL"], len(c["frag"]), fragment_seq=c["frag"],
                                config=kw.pop("config", PipelineConfig()), **kw)

    def test_junction_primers_span_the_junctions(self, construct):
        d = self._design(construct)
        jl, jr = d.amplicons[0], d.amplicons[1]
        JL, JR = construct["JL"], construct["JL"] + len(construct["frag"])
        assert (jl.role, jr.role) == ("junction-left", "junction-right")
        # forward primer wholly in the vector, reverse wholly in the fragment (left junction)
        assert jl.forward.end <= JL <= jl.reverse.start
        # forward primer wholly in the fragment, reverse wholly in the vector (right junction)
        assert jr.forward.end <= JR <= jr.reverse.start

    def test_every_product_is_reproduced_by_insilico_pcr(self, construct):
        d = self._design(construct)
        assert d.amplicons and not d.failed
        for a in d.amplicons:
            prods = predict_amplicons({"F": a.forward.sequence, "R": a.reverse.sequence},
                                      construct["seq"], circular=True)
            assert [(p["start"], p["size"]) for p in prods] == [(a.start, a.product_size)], a.name

    def test_empty_vector_gives_no_junction_product(self, construct):
        d = self._design(construct)
        empty = construct["v5"] + construct["v3"]           # re-circularised vector, no insert
        for a in d.amplicons[:2]:
            assert predict_amplicons({"F": a.forward.sequence, "R": a.reverse.sequence},
                                     empty, circular=True) == []

    def test_sizes_follow_distinct_ladder_ranges(self, construct):
        d = self._design(construct)
        for a, (lo, hi) in zip(d.amplicons, SIZE_LADDER):
            assert lo <= a.product_size <= hi, a.name
        sizes = [a.product_size for a in d.amplicons]
        assert len(set(sizes)) == len(sizes)

    def test_tm_in_configured_window_and_matches_pipeline_calc_tm(self, construct):
        cfg = PipelineConfig()
        d = self._design(construct, config=cfg)
        for a in d.amplicons:
            if a.level != "strict":
                continue
            for p in (a.forward, a.reverse):
                assert cfg.primer_tm_min - 0.15 <= p.tm <= cfg.primer_tm_max + 0.15
                assert p.tm == round(primer3.calc_tm(p.sequence), 1)

    def test_spaced_fallback_without_annotation(self, construct):
        d = self._design(construct)
        roles = [a.role for a in d.amplicons]
        assert roles == ["junction-left", "junction-right", "integrity", "integrity", "integrity"]
        assert any("no annotation" in n for n in d.notes)

    def test_n_spaced_is_configurable(self, construct):
        cfg = PipelineConfig(screening_n_spaced=1)
        d = self._design(construct, config=cfg)
        assert [a.role for a in d.amplicons].count("integrity") == 2    # 1 + 1 replacing the marker

    def test_integrity_amplicons_cover_distinct_regions(self, construct):
        d = self._design(construct)
        mids = [(a.start + a.end) / 2 for a in d.amplicons if a.role == "integrity"]
        assert mids == sorted(mids) and mids[-1] - mids[0] > 5000

    def test_inconsistent_coordinates_do_not_crash(self, construct):
        d = design_screening(construct["seq"], construct["JL"], 10 ** 7, config=PipelineConfig())
        assert d.amplicons == [] and d.notes

    def test_host_background_flags_planted_site(self, construct):
        d0 = self._design(construct)
        target = d0.amplicons[2].forward.sequence
        host = rnd(200_000, 9) + target + rnd(1000, 10)         # same primer sits in the host genome
        d = self._design(construct, host_seq=host)
        a = d.amplicons[2]
        assert a.forward.sequence != target                 # the host-contaminated primer was avoided
        assert not any("host genome" in w for w in a.warnings)

    def test_primer_with_second_site_on_construct_is_avoided(self):
        v5, v3 = rnd(3000, 21), rnd(2500, 22)
        frag = rnd(37000, 23, gc=True)
        cons = v5 + frag + v3
        base = design_screening(cons, len(v5), len(frag), config=PipelineConfig())
        first = base.amplicons[2]
        # Plant a perfect copy of that amplicon's forward primer far from every amplicon
        # window (so primer3's candidate list is unchanged): only the second-site check
        # can tell the two primers apart.
        pos = len(v5) + 15000
        assert all(not (a.start - 700 < pos < a.end + 700) for a in base.amplicons)
        plant = first.forward.sequence
        cons2 = cons[:pos] + plant + cons[pos + len(plant):]
        d2 = design_screening(cons2, len(v5), len(frag), config=PipelineConfig())
        a2 = d2.amplicons[2]
        assert a2.forward.sequence != first.forward.sequence
        assert not any("additional binding" in w for w in a2.warnings)


# ------------------------------------------------------------ marker gene
def _record(gene, strand, pad5=500, pad3=500, product="type II PKS ketosynthase alpha"):
    seq = rnd(pad5, 31) + (gene if strand == 1 else rc(gene)) + rnd(pad3, 32)
    loc = FeatureLocation(pad5, pad5 + len(gene), strand=strand)
    feats = [SeqFeature(loc, type="CDS", qualifiers={"product": [product], "locus_tag": ["KSA_01"]}),
             SeqFeature(FeatureLocation(10, 400, strand=1), type="CDS",
                        qualifiers={"product": ["hypothetical protein"], "locus_tag": ["HYP_01"]})]
    return SeqRecord(Seq(seq), id="bgc", features=feats)


class TestMarkerGene:
    GENE = "ATG" + rnd(1200, 41, gc=True) + "TGA"
    KW = ["t2pks", "ketosynthase", "chain length factor"]

    @pytest.mark.parametrize("strand", [1, -1])
    def test_found_on_either_strand_by_sequence(self, strand):
        rec = _record(self.GENE, strand)
        frag = str(rec.seq)
        m = find_marker_gene(rec, frag, self.KW)
        assert m and m["name"] == "KSA_01" and m["keyword"] == "ketosynthase"
        assert (m["end"] - m["start"]) == len(self.GENE)
        piece = frag[m["start"]:m["end"]]
        assert piece == (self.GENE if m["strand"] == "+" else rc(self.GENE))

    def test_located_in_a_fragment_that_is_shifted(self):
        rec = _record(self.GENE, 1)
        frag = rnd(777, 51) + str(rec.seq)                      # fragment carries extra flank
        m = find_marker_gene(rec, frag, self.KW)
        assert m["start"] == 777 + 500

    def test_no_match_or_no_features_returns_none(self):
        rec = _record(self.GENE, 1, product="hypothetical protein")
        assert find_marker_gene(rec, str(rec.seq), self.KW) is None
        assert find_marker_gene(None, "ACGT" * 100, self.KW) is None
        assert find_marker_gene(SeqRecord(Seq("ACGT" * 100)), "ACGT" * 100, self.KW) is None

    def test_keyword_priority_order(self):
        gene2 = "ATG" + rnd(900, 42) + "TAA"
        rec = _record(self.GENE, 1)
        seq = str(rec.seq) + rnd(300, 43) + gene2
        start2 = len(rec.seq) + 300
        rec2 = SeqRecord(Seq(seq), id="b", features=list(rec.features) + [
            SeqFeature(FeatureLocation(start2, start2 + len(gene2), strand=1), type="CDS",
                       qualifiers={"product": ["chain length factor"], "locus_tag": ["CLF_01"]})])
        assert find_marker_gene(rec2, seq, ["chain length factor", "ketosynthase"])["name"] == "CLF_01"
        assert find_marker_gene(rec2, seq, ["ketosynthase", "chain length factor"])["name"] == "KSA_01"

    def test_design_uses_marker_inside_the_gene(self):
        v5, v3 = rnd(3000, 61), rnd(2500, 62)
        frag = rnd(6000, 63, gc=True) + self.GENE + rnd(6000, 64, gc=True)
        rec = SeqRecord(Seq(frag), id="bgc", features=[SeqFeature(
            FeatureLocation(6000, 6000 + len(self.GENE), strand=1), type="CDS",
            qualifiers={"product": ["ketosynthase"], "locus_tag": ["KSA_01"]})])
        d = design_screening(v5 + frag + v3, len(v5), len(frag), fragment_seq=frag,
                             bgc_record=rec, config=PipelineConfig())
        mk = [a for a in d.amplicons if a.role == "marker-gene"]
        assert mk and d.marker["name"] == "KSA_01"
        gs, ge = len(v5) + 6000, len(v5) + 6000 + len(self.GENE)
        assert gs - 1 <= mk[0].forward.start and mk[0].reverse.end <= ge + 61
        assert [a.role for a in d.amplicons].count("integrity") == 2


# ------------------------------------------------------------ export
def test_csv_export(construct, tmp_path):
    d = design_screening(construct["seq"], construct["JL"], len(construct["frag"]),
                         fragment_seq=construct["frag"], config=PipelineConfig())
    out = export_screening(d, str(tmp_path / "s.csv"))
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 2 * len(d.amplicons)
    assert {"primer", "amplicon", "sequence_5to3", "tm_c", "product_bp", "partner"} <= set(rows[0])
    assert rows[0]["partner"] == rows[1]["primer"]
    assert all(set(r["sequence_5to3"]) <= set("ACGT") for r in rows)


# ------------------------------------------------------------ relative arm GC
class TestRelativeArmGc:
    def test_default_limits_unchanged_without_reference(self):
        assert gc_limits(PipelineConfig(), None) == (65.0, 75.0)

    def test_low_gc_genomes_keep_absolute_limits(self):
        for ref in (40, 50, 60):
            assert gc_limits(PipelineConfig(), ref) == (65.0, 75.0)

    def test_high_gc_fragment_raises_limits(self):
        ideal, warn = gc_limits(PipelineConfig(), 72)
        assert (ideal, warn) == (77.0, 82.0)

    def test_ideal_never_exceeds_warning(self):
        for ref in range(30, 90, 5):
            ideal, warn = gc_limits(PipelineConfig(), ref)
            assert ideal <= warn

    def test_72pct_arm_not_flagged_in_72pct_fragment_but_flagged_absolutely(self):
        arm = ("GC" * 7 + "GA" * 3 + "CG" * 3 + "AT" * 3 + "GC" * 3 + "GA" * 2 + "TT" * 2)[:50]
        gc = 100 * sum(b in "GC" for b in arm) / len(arm)
        assert 66 < gc < 77
        assert validate_gc(arm, PipelineConfig(), None)[0] is False
        assert validate_gc(arm, PipelineConfig(), 72.0)[0] is True

    def test_relative_can_be_disabled(self):
        cfg = PipelineConfig(relative_arm_gc=False)
        assert gc_limits(cfg, 72) == (65.0, 75.0)

    def test_extreme_arm_still_flagged_in_high_gc_fragment(self):
        arm = "GC" * 25
        ok, issues = validate_gc(arm, PipelineConfig(), 72.0)
        assert not ok and any("exceeds warning" in i for i in issues)

    def test_design_homology_arms_uses_fragment_gc(self, monkeypatch):
        from tar_crispr import homology_arms as ha
        seen = {}
        monkeypatch.setattr(ha, "find_valid_arm",
                            lambda frag, end, cfg, g, v, max_attempts=None, reference_gc=None:
                            seen.setdefault(end, reference_gc))
        frag = rnd(2000, 71, gc=True)
        ha.design_homology_arms(frag, frag, rnd(500, 72), PipelineConfig())
        assert seen["left"] == pytest.approx(100 * sum(b in "GC" for b in frag) / len(frag))
        seen.clear()
        ha.design_homology_arms(frag, frag, rnd(500, 72), PipelineConfig(relative_arm_gc=False))
        assert seen["left"] is None


# ------------------------------------------------------------ CLI end-to-end
class TestCli:
    def _run(self, tmp_path, *extra):
        import subprocess
        import sys
        from pathlib import Path
        td = Path(__file__).resolve().parent.parent / "test_data"
        out = tmp_path / "out"
        proc = subprocess.run(
            [sys.executable, "-m", "tar_crispr.cli", "run",
             "--bgc", str(td / "synthetic_bgc.fasta"), "--vector", str(td / "synthetic_vector.fasta"),
             "--genome", str(td / "synthetic_genome.fasta"), "--genbank", str(td / "synthetic_bgc.gbk"),
             "--output", str(out), "--no-blast", *extra], capture_output=True, text=True)
        return proc, out

    def test_report_and_csv_have_screening(self, tmp_path):
        proc, out = self._run(tmp_path)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        rep = (out / "report.md").read_text()
        assert "### 6.1 Colony-PCR Screening Panel" in rep
        assert "## 7. BGC Schematic Map" in rep and "## 8. Pipeline Configuration" in rep
        rows = list(csv.DictReader(open(out / "screening_primers.csv")))
        assert len(rows) >= 6

    def test_no_screening_flag(self, tmp_path):
        proc, out = self._run(tmp_path, "--no-screening")
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "Colony-PCR Screening Panel" not in (out / "report.md").read_text()
        assert not (out / "screening_primers.csv").exists()
