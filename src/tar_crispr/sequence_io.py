"""Step 1: Ingestion, validation, and cluster coordinate extraction."""
import glob
import os
import shutil
import subprocess
import tempfile
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


def locate_bgc_in_genome(bgc_record: SeqRecord,
                         flat_record: SeqRecord,
                         offsets: dict[str, tuple[int, int]],
                         evalue: float = 1e-10,
                         min_identity: float = 95.0,
                         min_coverage: float = 90.0,
                         ) -> tuple[str, int, int, float, list[str]]:
    """Locate *bgc_record* inside the flattened genome via BLAST.

    Used when the user supplies a genome (single- or multi-scaffold) and a
    BGC record but no explicit --start/--end: this finds where the BGC
    lives in the genome automatically, the same way we did manually with
    blastn + makeblastdb against scaffold2.fasta for the colinomycin case.

    Only ONE blastn database is built (against the whole flat genome),
    reused for the single query — this is cheap regardless of genome size,
    unlike the earlier per-candidate specificity bug.

    Returns
    -------
    (scaffold_id, local_start, local_end, identity_pct, warnings)
        Coordinates are 0-based, half-open, relative to the ORIGINAL
        scaffold (already translated back via translate_from_flat) — ready
        to pass straight into extract_cluster_bounds.

    Raises
    ------
    ValueError
        If no hit clears min_identity/min_coverage, or every hit spans a
        scaffold separator (spurious cross-scaffold alignment).
    """
    warnings: list[str] = []
    bgc_seq = str(bgc_record.seq).upper()
    bgc_len = len(bgc_seq)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qfh, \
         tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as dfh:
        qfh.write(f">query\n{bgc_seq}\n")
        qpath = qfh.name
        dfh.write(f">{flat_record.id}\n{str(flat_record.seq).upper()}\n")
        dpath = dfh.name

    db_dir = os.path.dirname(dpath)
    db_name = os.path.splitext(os.path.basename(dpath))[0]
    full_db = os.path.join(db_dir, db_name)

    try:
        subprocess.run(
            ["makeblastdb", "-in", dpath, "-dbtype", "nucl", "-out", full_db],
            check=True, capture_output=True, timeout=120,
        )
        result = subprocess.run(
            ["blastn", "-query", qpath, "-db", full_db,
             "-evalue", str(evalue), "-max_hsps", "10", "-dust", "no",
             "-outfmt", "6 sseqid pident length qlen sstart send bitscore"],
            check=True, capture_output=True, text=True, timeout=120,
        )
    finally:
        # BLAST+ database file extensions vary by version (v4 vs v5 index
        # format: .nhr/.nin/.nsq vs .ndb/.njs/.not/.ntf/.nto/...). Glob for
        # anything matching the db prefix instead of a fixed extension list,
        # so this doesn't silently leak temp files on a version bump.
        for p in glob.glob(full_db + "*"):
            if os.path.exists(p):
                os.unlink(p)
        for p in (dpath, qpath):
            if os.path.exists(p):
                os.unlink(p)

    hits = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        sseqid, pident, length, qlen, sstart, send, bitscore = line.split("\t")
        hits.append({
            "pident": float(pident), "length": int(length), "qlen": int(qlen),
            "sstart": int(sstart), "send": int(send), "bitscore": float(bitscore),
        })
    if not hits:
        raise ValueError(
            f"BGC not found in genome (no BLAST hits at evalue<={evalue})."
        )

    hits.sort(key=lambda h: h["bitscore"], reverse=True)

    best_score = hits[0]["bitscore"]
    close_hits = [h for h in hits if h["bitscore"] >= 0.5 * best_score]
    if len(close_hits) > 1:
        warnings.append(
            f"{len(close_hits)} BLAST hits with similar score found "
            f"(best={best_score:.0f}); picked the top-scoring one — verify "
            f"this is the correct locus if the genome has repeats."
        )

    for h in hits:
        s, e = h["sstart"], h["send"]
        if s > e:
            s, e = e, s
            warnings.append("BGC hit on the reverse strand of the genome.")
        flat_start, flat_end = s - 1, e  # BLAST is 1-based inclusive

        try:
            scaffold_id, local_start, local_end = translate_from_flat(
                flat_start, flat_end, offsets
            )
        except ValueError:
            continue  # spurious cross-scaffold alignment, try next hit

        coverage = 100.0 * h["length"] / bgc_len if bgc_len else 0.0
        if h["pident"] < min_identity:
            warnings.append(
                f"Identity {h['pident']:.1f}% is below the {min_identity}% "
                f"threshold — this may not be the exact source locus."
            )
        if coverage < min_coverage:
            warnings.append(
                f"Coverage {coverage:.1f}% is below the {min_coverage}% "
                f"threshold — the BGC may be split across scaffolds or "
                f"only partially present in this genome."
            )
        if h["pident"] < min_identity or coverage < min_coverage:
            continue  # keep looking for a better hit before giving up

        return scaffold_id, local_start, local_end, h["pident"], warnings

    raise ValueError(
        f"No BLAST hit for the BGC cleared identity>={min_identity}% and "
        f"coverage>={min_coverage}% without spanning a scaffold separator. "
        f"Provide --start/--end manually."
    )


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
