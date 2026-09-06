"""Step 1: Ingestion, validation, and cluster coordinate extraction."""
import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord

from tar_crispr.config import PipelineConfig, ClusterInfo


def read_sequence(path: str) -> SeqRecord:
    """Read a sequence from FASTA or GenBank file, auto-detecting format.

    Parameters
    ----------
    path : str
        Path to the sequence file (.fasta, .fa, .fna, .gb, .gbk, .genbank).

    Returns
    -------
    SeqRecord
        Biopython SeqRecord with the sequence and metadata.
    """
    fmt = _detect_format(path)
    records = list(SeqIO.parse(path, fmt))
    if len(records) == 0:
        raise ValueError(f"No sequences found in {path}")
    if len(records) > 1:
        return SeqRecord(Seq("".join(str(r.seq) for r in records)),
                         id="concatenated",
                         description=";".join(r.id for r in records),
                         annotations={"molecule_type": "DNA"})
    rec = records[0]
    return rec


def _detect_format(path: str) -> str:
    """Detect file format from extension or content."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".fasta", ".fa", ".fna"):
        return "fasta"
    if ext in (".gb", ".gbk", ".genbank"):
        return "genbank"
    # Sniff content
    with open(path, "r") as fh:
        first_line = fh.readline().strip()
    if first_line.startswith("LOCUS"):
        return "genbank"
    if first_line.startswith(">"):
        return "fasta"
    raise ValueError(f"Cannot determine format of {path}")


def validate_sequence(seq: SeqRecord) -> dict:
    """Validate sequence integrity and return stats.

    Parameters
    ----------
    seq : SeqRecord
        The sequence to validate.

    Returns
    -------
    dict
        Dictionary with 'length', 'gc_percent', 'invalid_chars', 'valid'.
    """
    seq_str = str(seq.seq).upper()
    valid_nucleotides = set("ACGTN")
    invalid = [c for c in seq_str if c not in valid_nucleotides]
    stats = {
        "length": len(seq_str),
        "gc_percent": round(gc_fraction(seq_str) * 100, 2) if len(seq_str) > 0 else 0.0,
        "invalid_chars": sorted(set(invalid)),
        "valid": len(invalid) == 0 and len(seq_str) > 0,
        "seqrecord_id": seq.id,
        "description": seq.description,
    }
    return stats


def extract_cluster_bounds(seqrecord: SeqRecord,
                           start: Optional[int] = None,
                           end: Optional[int] = None) -> ClusterInfo:
    """Identify cluster start/end coordinates.

    If start and end are provided, use them directly. Otherwise, parse from
    GenBank features (looks for 'cluster' feature type or 'CDS' span).

    Parameters
    ----------
    seqrecord : SeqRecord
        The full genomic sequence record.
    start : int, optional
        1-based start coordinate of the BGC.
    end : int, optional
        1-based end coordinate of the BGC.

    Returns
    -------
    ClusterInfo
        Dataclass with cluster coordinates.
    """
    if start is not None and end is not None:
        return ClusterInfo(start=start, end=end)

    # Try parsing GenBank features
    if seqrecord.features:
        # Look for 'cluster' feature
        for f in seqrecord.features:
            if f.type.lower() == "cluster":
                s = int(f.location.start)
                e = int(f.location.end)
                return ClusterInfo(start=s, end=e,
                                   name=f.qualifiers.get("label", ["BGC"])[0])
        # Fallback: use CDS span
        cds = [f for f in seqrecord.features if f.type.lower() == "cds"]
        if cds:
            s = min(int(f.location.start) for f in cds)
            e = max(int(f.location.end) for f in cds)
            return ClusterInfo(start=s, end=e, name="CDS_span")

    raise ValueError(
        "Cluster bounds not provided and could not be parsed from GenBank features. "
        "Please provide --start and --end."
    )


def get_flanking_sequence(full_seq: SeqRecord, cluster: ClusterInfo,
                          window: int = 500) -> tuple[str, str, str]:
    """Extract the upstream flank, the cluster, and downstream flank.

    Parameters
    ----------
    full_seq : SeqRecord
        Full genomic sequence.
    cluster : ClusterInfo
        Cluster coordinates.
    window : int
        Size of flanking window in bp (default 500).

    Returns
    -------
    tuple of (upstream_seq, cluster_seq, downstream_seq)
    """
    seq_str = str(full_seq.seq).upper()
    c_start = max(0, cluster.start - window)
    c_end = cluster.end
    upstream = seq_str[c_start:cluster.start]
    cluster_seq = seq_str[cluster.start:cluster.end]
    downstream_end = min(len(seq_str), cluster.end + window)
    downstream = seq_str[cluster.end:downstream_end]
    return upstream, cluster_seq, downstream


def check_blast_available() -> bool:
    """Check if blastn is available on the system."""
    return shutil.which("blastn") is not None


def check_rnafold_available() -> bool:
    """Check if RNAfold is available on the system."""
    return shutil.which("RNAfold") is not None


def check_external_tools(config: PipelineConfig) -> PipelineConfig:
    """Detect external tools and update config. Returns updated config."""
    config.blast_available = check_blast_available()
    config.rnafold_available = check_rnafold_available()
    return config
