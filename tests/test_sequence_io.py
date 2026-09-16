"""Tests for sequence_io.py (Step 1)."""
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tar_crispr.sequence_io import (
    read_sequence,
    read_genome,
    translate_to_flat,
    translate_from_flat,
    locate_bgc_in_genome,
    validate_sequence,
    extract_cluster_bounds,
    get_flanking_sequence,
    check_blast_available,
    check_external_tools,
)
from tar_crispr.config import PipelineConfig
import shutil


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


class TestReadGenome:
    def test_single_record(self, tmp_path):
        rec = SeqRecord(Seq("ACGT" * 50), id="only_scaffold", description="")
        path = tmp_path / "single.fasta"
        with open(path, "w") as fh:
            SeqIO.write(rec, fh, "fasta")

        flat, offsets = read_genome(str(path))
        assert str(flat.seq) == "ACGT" * 50
        assert offsets == {"only_scaffold": (0, 200)}

    def test_multi_record_has_separators(self, tmp_path):
        recs = [
            SeqRecord(Seq("A" * 100), id="scaf1", description=""),
            SeqRecord(Seq("C" * 150), id="scaf2", description=""),
            SeqRecord(Seq("G" * 80), id="scaf3", description=""),
        ]
        path = tmp_path / "multi.fasta"
        with open(path, "w") as fh:
            SeqIO.write(recs, fh, "fasta")

        flat, offsets = read_genome(str(path), sep_len=60)

        assert offsets["scaf1"] == (0, 100)
        assert offsets["scaf2"] == (160, 310)
        assert offsets["scaf3"] == (370, 450)
        assert len(flat.seq) == 100 + 60 + 150 + 60 + 80
        assert str(flat.seq[100:160]) == "N" * 60

    def test_duplicate_ids_get_suffixed(self, tmp_path):
        recs = [
            SeqRecord(Seq("A" * 10), id="dup", description=""),
            SeqRecord(Seq("C" * 10), id="dup", description=""),
        ]
        path = tmp_path / "dup.fasta"
        with open(path, "w") as fh:
            SeqIO.write(recs, fh, "fasta")

        flat, offsets = read_genome(str(path))
        assert len(offsets) == 2
        assert "dup" in offsets
        assert any(k != "dup" for k in offsets)


class TestTranslateCoordinates:
    def test_round_trip(self, tmp_path):
        recs = [
            SeqRecord(Seq("A" * 500), id="s1", description=""),
            SeqRecord(Seq("C" * 700), id="s2", description=""),
        ]
        path = tmp_path / "genome.fasta"
        with open(path, "w") as fh:
            SeqIO.write(recs, fh, "fasta")
        flat, offsets = read_genome(str(path))

        flat_start, flat_end = translate_to_flat("s2", 50, 200, offsets)
        scaffold_id, local_start, local_end = translate_from_flat(
            flat_start, flat_end, offsets
        )
        assert scaffold_id == "s2"
        assert (local_start, local_end) == (50, 200)

    def test_from_flat_raises_across_separator(self, tmp_path):
        recs = [
            SeqRecord(Seq("A" * 100), id="s1", description=""),
            SeqRecord(Seq("C" * 100), id="s2", description=""),
        ]
        path = tmp_path / "genome.fasta"
        with open(path, "w") as fh:
            SeqIO.write(recs, fh, "fasta")
        flat, offsets = read_genome(str(path), sep_len=60)

        with pytest.raises(ValueError, match="separator"):
            translate_from_flat(90, 170, offsets)


@pytest.mark.skipif(not shutil.which("blastn"), reason="blastn not installed")
class TestLocateBgcInGenome:
    def _make_genome_with_bgc(self, tmp_path, bgc_seq, scaffold_index=1,
                              n_scaffolds=3, flank=300):
        import random
        random.seed(7)
        parts = []
        bgc_scaffold_id = None
        for i in range(n_scaffolds):
            if i == scaffold_index:
                filler = "".join(random.choices("ACGT", k=flank))
                seq = filler + bgc_seq + filler
                sid = f"scaf{i}_with_bgc"
                bgc_scaffold_id = sid
            else:
                seq = "".join(random.choices("ACGT", k=flank * 2))
                sid = f"scaf{i}"
            parts.append(SeqRecord(Seq(seq), id=sid, description=""))

        genome_path = tmp_path / "genome.fasta"
        with open(genome_path, "w") as fh:
            SeqIO.write(parts, fh, "fasta")
        return str(genome_path), bgc_scaffold_id, flank

    def test_finds_bgc_in_correct_scaffold(self, tmp_path):
        import random
        random.seed(42)
        # Non-repetitive synthetic BGC (tandem repeats confuse BLAST's HSP
        # extension — real BGCs aren't purely repetitive either).
        bgc_seq = "".join(random.choices("ACGT", k=600))
        genome_path, expected_scaffold, flank = self._make_genome_with_bgc(
            tmp_path, bgc_seq, scaffold_index=1, n_scaffolds=4
        )
        bgc_path = tmp_path / "bgc.fasta"
        with open(bgc_path, "w") as fh:
            SeqIO.write(SeqRecord(Seq(bgc_seq), id="query_bgc", description=""),
                       fh, "fasta")

        bgc = read_sequence(str(bgc_path))
        flat, offsets = read_genome(genome_path)

        scaffold_id, start, end, identity, warnings = locate_bgc_in_genome(
            bgc, flat, offsets
        )

        assert scaffold_id == expected_scaffold
        assert start == flank
        assert end == flank + len(bgc_seq)
        assert identity == pytest.approx(100.0, abs=0.01)
        assert warnings == []

    def test_raises_when_not_present(self, tmp_path):
        genome_path, _, _ = self._make_genome_with_bgc(
            tmp_path, "ACGT" * 100, scaffold_index=0, n_scaffolds=2
        )
        unrelated_bgc = "".join(["TTTTAAAACCCCGGGG"] * 20)
        bgc_path = tmp_path / "bgc.fasta"
        with open(bgc_path, "w") as fh:
            SeqIO.write(SeqRecord(Seq(unrelated_bgc), id="query", description=""),
                       fh, "fasta")

        bgc = read_sequence(str(bgc_path))
        flat, offsets = read_genome(genome_path)

        with pytest.raises(ValueError):
            locate_bgc_in_genome(bgc, flat, offsets)

    def test_cleans_up_all_temp_blast_files(self, tmp_path):
        # Regression: the old cleanup used a fixed extension list
        # (.nhr/.nin/.nsq/...) from BLAST+'s v4 db format, but BLAST+ 2.10+
        # uses v5 (.ndb/.njs/.not/.ntf/.nto/...) and leaked 5 files per call
        # into /tmp. Cleanup now globs the db prefix instead of guessing
        # extensions.
        import glob
        bgc_seq = "".join(["GCGTACGATCGATCGATGCGCGATCGATGCATGCATGC"] * 5)
        genome_path, _, _ = self._make_genome_with_bgc(
            tmp_path, bgc_seq, scaffold_index=0, n_scaffolds=2
        )
        bgc_path = tmp_path / "bgc.fasta"
        with open(bgc_path, "w") as fh:
            SeqIO.write(SeqRecord(Seq(bgc_seq), id="query", description=""),
                       fh, "fasta")
        bgc = read_sequence(str(bgc_path))
        flat, offsets = read_genome(genome_path)

        before = set(glob.glob("/tmp/tmp*"))
        locate_bgc_in_genome(bgc, flat, offsets)
        after = set(glob.glob("/tmp/tmp*"))

        leftover = after - before
        assert not leftover, (
            f"locate_bgc_in_genome leaked {len(leftover)} temp file(s): "
            f"{sorted(leftover)}"
        )
