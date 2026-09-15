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


def read_genome(path: str, sep_len: int = 60
                ) -> tuple[SeqRecord, dict[str, tuple[int, int]]]:
    """Read a (possibly multi-record) genome file.

    Returns a single flattened SeqRecord that concatenates every record in
    the file, separated by ``sep_len`` Ns, plus a map from scaffold ID to
    the (start, end) half-open interval it occupies in the flat sequence.

    The N separator prevents 20 nt sgRNAs (or 50-60 bp homology arms) from
    spuriously matching across scaffold boundaries when the flat genome is
    scanned for specificity or uniqueness.

    Parameters
    ----------
    path : str
        Path to a FASTA or GenBank file (single- or multi-record).
    sep_len : int
        Length of the N separator inserted between records (default 60).

    Returns
    -------
    flat : SeqRecord
        Concatenated genome with separators, id="flattened_genome".
    offsets : dict
        Mapping ``scaffold_id -> (start, end)`` in the flat sequence.
    """
    fmt = _detect_format(path)
    records = list(SeqIO.parse(path, fmt))
    if not records:
        raise ValueError(f"No sequences found in {path}")

    sep = "N" * sep_len
    parts: list[str] = []
    offsets: dict[str, tuple[int, int]] = {}
    cursor = 0

    for i, rec in enumerate(records):
        sid = rec.id or f"record_{i}"
        # Handle duplicate IDs by appending a suffix
        if sid in offsets:
            sid = f"{sid}_{i}"
        seq_str = str(rec.seq).upper()
        seq_len = len(seq_str)
        offsets[sid] = (cursor, cursor + seq_len)
        parts.append(seq_str)
        cursor += seq_len
        # Do NOT add separator after the last record
        if i < len(records) - 1:
            parts.append(sep)
            cursor += sep_len

    flat_seq = "".join(parts)
    flat_rec = SeqRecord(
        Seq(flat_seq),
        id="flattened_genome",
        description=f"{len(records)} scaffolds; separator={sep_len}N",
        annotations={"molecule_type": "DNA"},
    )
    return flat_rec, offsets


def translate_to_flat(scaffold_id: str, local_start: int, local_end: int,
                      offsets: dict[str, tuple[int, int]]) -> tuple[int, int]:
    """Translate (start, end) relative to a scaffold into coordinates in
    the flattened genome.

    Parameters
    ----------
    scaffold_id : str
        Scaffold identifier as stored in *offsets*.
    local_start, local_end : int
        Half-open coordinates relative to the scaffold.
    offsets : dict
        Map returned by :func:`read_genome`.

    Returns
    -------
    tuple of (flat_start, flat_end)
    """
    if scaffold_id not in offsets:
        raise KeyError(
            f"Scaffold {scaffold_id!r} not found in genome. "
            f"Available: {sorted(offsets)}"
        )
    base, _ = offsets[scaffold_id]
    return base + local_start, base + local_end


def translate_from_flat(flat_start: int, flat_end: int,
                        offsets: dict[str, tuple[int, int]]
                        ) -> tuple[str, int, int]:
    """Map a coordinate in the flattened genome back to (scaffold, local_start, local_end).

    Parameters
    ----------
    flat_start, flat_end : int
        Half-open coordinates in the flattened genome.
    offsets : dict
        Map returned by :func:`read_genome`.

    Returns
    -------
    (scaffold_id, local_start, local_end)
    """
    for sid, (s, e) in offsets.items():
        if s <= flat_start and flat_end <= e:
            return sid, flat_start - s, flat_end - s
    raise ValueError(
        f"Coordinates {flat_start}..{flat_end} span a separator or fall "
        f"outside every scaffold."
    )


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
                           end: Optional[int] = None,
                           gene_kinds: Optional[list[str]] = None) -> ClusterInfo:
    """Identify cluster start/end coordinates.

    Resolution order:

    1. If ``start`` and ``end`` are given, use them directly.
    2. If the record has an antiSMASH ``region`` feature, use it:

       - with ``gene_kinds`` set (e.g. ``["biosynthetic",
         "biosynthetic-additional"]``), the bounds are computed from the
         span of the CDS inside the region whose ``gene_kind`` qualifier
         is in the list;
       - with ``gene_kinds`` unset or ``["all"]``, the bounds are the
         start/end of the region feature itself.
    3. Otherwise, look for a generic ``cluster`` feature.
    4. Otherwise, use the span of all ``CDS`` features.

    Parameters
    ----------
    seqrecord : SeqRecord
        The genomic sequence record (from FASTA or GenBank).
    start, end : int, optional
        Explicit 0-based half-open coordinates. Take precedence over
        any annotation.
    gene_kinds : list of str, optional
        Restrict antiSMASH region bounds to CDS matching these
        ``gene_kind`` values. Use ``["all"]`` to force the region span.

    Returns
    -------
    ClusterInfo
        Dataclass with cluster coordinates.
    """
    if start is not None and end is not None:
        return ClusterInfo(start=start, end=end)

    # Try antiSMASH region first
    if seqrecord.features:
        region = _find_antismash_region(seqrecord)
        if region is not None:
            return _bounds_from_antismash_region(region, seqrecord, gene_kinds)

        # Generic 'cluster' feature
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


def _find_antismash_region(seqrecord: SeqRecord):
    """Return the first antiSMASH 'region' feature, or None."""
    for f in seqrecord.features:
        if f.type.lower() == "region":
            return f
    return None


def _bounds_from_antismash_region(region, seqrecord: SeqRecord,
                                  gene_kinds: Optional[list[str]]) -> ClusterInfo:
    """Compute cluster bounds from an antiSMASH region feature.

    If ``gene_kinds`` is None, empty, or contains only "all", the region
    span is used directly. Otherwise the CDS inside the region are
    filtered by their ``gene_kind`` qualifier, and the bounds are the
    min(start) / max(end) of the filtered CDS.
    """
    region_start = int(region.location.start)
    region_end = int(region.location.end)
    region_number = region.qualifiers.get("region_number", ["?"])[0]
    product = region.qualifiers.get("product", [""])[0]
    name = f"antiSMASH_region_{region_number}"
    description = product

    kinds = None
    if gene_kinds:
        # Accept both "--gene-kinds a,b" (single string with commas) and
        # "--gene-kinds a --gene-kinds b" (list of strings).
        flat: list[str] = []
        for k in gene_kinds:
            flat.extend(k.split(","))
        normalized = [x.strip().lower() for x in flat if x.strip()]
        if normalized and normalized != ["all"]:
            kinds = set(normalized)

    if kinds is None:
        return ClusterInfo(
            start=region_start, end=region_end,
            name=name, description=description,
        )

    # Filter CDS inside the region by gene_kind
    selected = []
    for f in seqrecord.features:
        if f.type.lower() != "cds":
            continue
        s = int(f.location.start)
        e = int(f.location.end)
        if s < region_start or e > region_end:
            continue
        kind = f.qualifiers.get("gene_kind", [""])[0].lower()
        if kind in kinds:
            selected.append((s, e))

    if not selected:
        raise ValueError(
            f"No CDS inside antiSMASH region with gene_kind in {sorted(kinds)}. "
            f"Available kinds: {_available_gene_kinds(seqrecord)}"
        )

    s = min(c[0] for c in selected)
    e = max(c[1] for c in selected)
    return ClusterInfo(
        start=s, end=e, name=name,
        description=f"{description} (filtered by {sorted(kinds)})",
    )


def _available_gene_kinds(seqrecord: SeqRecord) -> list[str]:
    """Return the sorted list of gene_kind values present on CDS features."""
    kinds = set()
    for f in seqrecord.features:
        if f.type.lower() != "cds":
            continue
        for k in f.qualifiers.get("gene_kind", []):
            kinds.add(k.lower())
    return sorted(kinds)


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
