"""Tests for sequence_io.py (Step 1)."""
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tar_crispr.sequence_io import (
    read_sequence,
    validate_sequence,
    extract_cluster_bounds,
    get_flanking_sequence,
    check_blast_available,
    check_external_tools,
)
from tar_crispr.config import PipelineConfig


class TestReadSequence:
    def test_read_fasta(self, tmp_fasta):
        path = tmp_fasta("test.fasta", "ACGTACGTACGT")
        rec = read_sequence(path)
        assert str(rec.seq) == "ACGTACGTACGT"

    def test_read_genbank(self, tmp_gbk):
        path = tmp_gbk("test.gb", "ACGTACGTACGTACGTACGT", 5, 15)
        rec = read_sequence(path)
        assert str(rec.seq) == "ACGTACGTACGTACGTACGT"
        assert len(rec.features) > 0

    def test_read_multifasta(self, tmp_path):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        recs = [
            SeqRecord(Seq("AAAA"), id="s1", description=""),
            SeqRecord(Seq("CCCC"), id="s2", description=""),
        ]
        path = str(tmp_path / "multi.fasta")
        with open(path, "w") as fh:
            SeqIO.write(recs, fh, "fasta")
        rec = read_sequence(path)
        assert str(rec.seq) == "AAAACCCC"

    def test_read_nonexistent(self):
        with pytest.raises((FileNotFoundError, Exception)):
            read_sequence("/nonexistent/path/file.fasta")


class TestValidateSequence:
    def test_valid_sequence(self):
        rec = SeqRecord(Seq("ACGTACGTACGT"), id="test")
        stats = validate_sequence(rec)
        assert stats["valid"] is True
        assert stats["length"] == 12
        assert stats["gc_percent"] > 0

    def test_invalid_characters(self):
        rec = SeqRecord(Seq("ACGTXACGT"), id="test")
        stats = validate_sequence(rec)
        assert stats["valid"] is False
        assert "X" in stats["invalid_chars"]

    def test_gc_content(self):
        rec = SeqRecord(Seq("GCGCGCGCGCGC"), id="test")
        stats = validate_sequence(rec)
        assert stats["gc_percent"] == pytest.approx(100.0, abs=0.01)

    def test_empty_sequence(self):
        rec = SeqRecord(Seq(""), id="test")
        stats = validate_sequence(rec)
        assert stats["valid"] is False
        assert stats["length"] == 0

    def test_with_n_characters(self):
        rec = SeqRecord(Seq("ACGTNNNNACGT"), id="test")
        stats = validate_sequence(rec)
        assert stats["valid"] is True  # N is allowed


class TestExtractClusterBounds:
    def test_with_explicit_coords(self):
        rec = SeqRecord(Seq("A" * 1000), id="test")
        rec.features = []
        cluster = extract_cluster_bounds(rec, start=100, end=500)
        assert cluster.start == 100
        assert cluster.end == 500

    def test_from_genbank_cluster_feature(self, tmp_gbk):
        path = tmp_gbk("test.gb", "A" * 100, 20, 80)
        rec = read_sequence(path)
        cluster = extract_cluster_bounds(rec)
        assert cluster.start == 20
        assert cluster.end == 80
        assert cluster.name == "myBGC"

    def test_no_coords_no_features(self):
        rec = SeqRecord(Seq("A" * 100), id="test")
        rec.features = []
        with pytest.raises(ValueError):
            extract_cluster_bounds(rec)


class TestGetFlankingSequence:
    def test_flanking_sequence(self):
        seq_str = "A" * 500 + "B" * 1000 + "C" * 500
        rec = SeqRecord(Seq(seq_str), id="test")
        cluster = extract_cluster_bounds(rec, start=500, end=1500)
        upstream, cluster_seq, downstream = get_flanking_sequence(rec, cluster, window=500)
        assert upstream == "A" * 500
        assert cluster_seq == "B" * 1000
        assert downstream == "C" * 500

    def test_flanking_short_sequence(self):
        seq_str = "AAAA" + "BBBB" + "CCCC"
        rec = SeqRecord(Seq(seq_str), id="test")
        cluster = extract_cluster_bounds(rec, start=4, end=8)
        upstream, cluster_seq, downstream = get_flanking_sequence(rec, cluster, window=10)
        assert upstream == "AAAA"
        assert cluster_seq == "BBBB"
        assert downstream == "CCCC"


class TestExternalTools:
    def test_check_blast_available(self):
        result = check_blast_available()
        assert isinstance(result, bool)

    def test_check_external_tools(self):
        config = PipelineConfig()
        config = check_external_tools(config)
        assert isinstance(config.blast_available, bool)
        assert isinstance(config.rnafold_available, bool)
