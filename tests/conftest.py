"""Pytest fixtures with synthetic sequences for testing."""
import os
import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def project_root():
    return Path(__file__).parent.parent


@pytest.fixture
def tmp_fasta(tmp_path):
    """Create a temporary FASTA file."""
    def _make(name, sequence, seq_id="test"):
        path = tmp_path / name
        rec = SeqRecord(Seq(sequence.upper()), id=seq_id, description="")
        with open(path, "w") as fh:
            SeqIO.write(rec, fh, "fasta")
        return str(path)
    return _make


@pytest.fixture
def tmp_gbk(tmp_path):
    """Create a temporary GenBank file with a cluster feature."""
    def _make(name, sequence, cluster_start, cluster_end, seq_id="test"):
        path = tmp_path / name
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        rec = SeqRecord(Seq(sequence.upper()), id=seq_id, description="")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["genbank_division"] = "BCT"
        rec.annotations["date"] = "01-JAN-2026"
        # Cluster feature
        loc = FeatureLocation(cluster_start, cluster_end, strand=1)
        feat = SeqFeature(loc, type="cluster", qualifiers={"label": ["myBGC"]})
        rec.features.append(feat)
        with open(path, "w") as fh:
            SeqIO.write(rec, fh, "genbank")
        return str(path)
    return _make


@pytest.fixture
def synthetic_bgc_seq():
    """A ~1000bp synthetic BGC sequence with gene-like structure."""
    # Create a sequence with a realistic pattern: high GC, gene-like structure
    # Contains PAM sites (NGG) on both strands for testing
    parts = []
    # Promoter-like region with high GC
    parts.append("GCGCGCGGCGCGGCGCGGCGCGGCGCGGCGCGGCGCGGCGCGGCGCGGCGCGGCGCGGCGC")
    parts.append("ATGGCTGACGACGTCGACGTCGACGTCGACGTCGACGTCGACGTCGACGTCGACGTCGAC")
    parts.append("GGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")
    parts.append("CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA")
    parts.append("GGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")
    parts.append("CCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGG")
    parts.append("ATGCCGCTGACCGGATCCTGACCGGATCCGACCGGATCCGACCGGATCCGACCGGATCCGA")
    parts.append("TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG")
    parts.append("GGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCC")
    parts.append("AATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAA")
    parts.append("CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCG")
    parts.append("GGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")
    full = "".join(parts)
    # Pad to approximately 960bp
    padding = "AT" * 40
    full = full[:960] + padding
    return full[:1000] if len(full) > 1000 else full.ljust(1000, "N")


@pytest.fixture
def synthetic_genome_seq(synthetic_bgc_seq):
    """Genome = 1000bp upstream + BGC + 1000bp downstream."""
    upstream = "A" * 1000
    downstream = "T" * 1000
    return upstream + synthetic_bgc_seq + downstream, 1000, 1000 + len(synthetic_bgc_seq)


@pytest.fixture
def synthetic_vector_seq():
    """A ~600bp synthetic capture vector with a unique restriction site."""
    # Multiple cloning site with EcoRI site
    backbone = "G" * 300
    mcs = "GAATTCGCGGCCGCTCTAGAGGATCCACTAGTCGAGCCATGG"
    backbone2 = "C" * 296
    return backbone + mcs + backbone2


@pytest.fixture
def synthetic_full_test(tmp_path, synthetic_bgc_seq, synthetic_genome_seq, synthetic_vector_seq):
    """Full synthetic test data: writes all files to tmp_path."""
    genome_seq, bgc_start, bgc_end = synthetic_genome_seq
    files = {}

    bgc_rec = SeqRecord(Seq(genome_seq), id="synthetic_genome", description="")
    bgc_path = str(tmp_path / "synthetic_bgc.fasta")
    with open(bgc_path, "w") as fh:
        SeqIO.write(bgc_rec, fh, "fasta")
    files["bgc_fasta"] = bgc_path

    genome_rec = SeqRecord(Seq(genome_seq), id="synthetic_genome", description="")
    genome_path = str(tmp_path / "synthetic_genome.fasta")
    with open(genome_path, "w") as fh:
        SeqIO.write(genome_rec, fh, "fasta")
    files["genome_fasta"] = genome_path

    vec_rec = SeqRecord(Seq(synthetic_vector_seq), id="pCatch_vector", description="")
    vec_path = str(tmp_path / "synthetic_vector.fasta")
    with open(vec_path, "w") as fh:
        SeqIO.write(vec_rec, fh, "fasta")
    files["vector_fasta"] = vec_path

    files["bgc_start"] = bgc_start
    files["bgc_end"] = bgc_end
    return files


@pytest.fixture
def clean_env():
    """Set up clean environment for tests."""
    old_env = dict(os.environ)
    yield
    os.environ.clear()
    os.environ.update(old_env)
