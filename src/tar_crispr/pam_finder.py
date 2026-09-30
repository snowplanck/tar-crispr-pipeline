"""Step 2: PAM detection, sgRNA candidate identification, scoring, and ranking.

For SpCas9 (PAM = NGG):
  - Forward-strand sgRNA: protospacer(20nt) is 5' of an NGG on the forward strand.
    Cas9 cuts 3 bp upstream of the PAM (blunt cut), i.e. at position PAM_start - 3.
  - Reverse-strand sgRNA: PAM (NGG on reverse complement = CCN on forward)
    is followed on the forward strand by the reverse-complement protospacer.
    Cas9 cuts 3 bp upstream of the PAM on the reverse strand, which maps to
    position (len(seq) - (PAM_rc_start + protospacer_length + 3)) on the forward strand.
"""
import os
import re
import shutil
import subprocess
import tempfile
import dataclasses
from typing import Optional
from Bio.SeqUtils import gc_fraction
from tar_crispr.config import PipelineConfig, PAMCandidate
from tar_crispr.cut_specificity import build_site_index, find_cut_sites


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence."""
    comp_table = str.maketrans("ACGTNacgtn", "TGCANTGCAN")
    return seq.translate(comp_table)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return complement(seq)[::-1]


def _build_pam_regex(pam: str) -> str:
    """Convert a PAM IUPAC pattern (e.g. 'NGG') into a regex string."""
    iupac = {
        "N": "[ACGT]", "R": "[AG]", "Y": "[CT]", "M": "[AC]",
        "K": "[GT]", "S": "[CG]", "W": "[AT]", "H": "[ACT]",
        "B": "[CGT]", "V": "[ACG]", "D": "[AGT]",
    }
    return "".join(iupac.get(ch, ch) for ch in pam.upper())


def _check_polyt(sequence: str, threshold: int = 4) -> bool:
    """Flag sequences containing a run of *threshold* or more consecutive T's."""
    return "T" * threshold in sequence.upper()


def calculate_gc(seq: str) -> float:
    """Return GC percentage (0-100 scale)."""
    if len(seq) == 0:
        return 0.0
    return round(gc_fraction(seq) * 100, 2)


def find_pam_sites(sequence: str,
                   pam: str = "NGG",
                   protospacer_length: int = 20) -> list[PAMCandidate]:
    """Scan *sequence* on both strands for PAM sites and extract sgRNA candidates.

    Parameters
    ----------
    sequence : str
        Target DNA sequence.
    pam : str
        PAM IUPAC pattern (default ``"NGG"`` for SpCas9).
    protospacer_length : int
        Length of the guide RNA (default ``20``).

    Returns
    -------
    list[PAMCandidate]
        All candidate sgRNAs found on both strands.
    """
    seq = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")
    candidates: list[PAMCandidate] = []
    pam_regex = _build_pam_regex(pam)
    half = protospacer_length

    # ---- Forward strand ----
    # Pattern: PROTOSPACER(20nt) + PAM(NGG)
    # We search for PAM, then grab the 20 nt immediately upstream.
    for m in re.finditer(pam_regex, seq):
        pam_start = m.start()
        proto_start = pam_start - half
        if proto_start < 0:
            continue
        protospacer = seq[proto_start:pam_start]
        # Cas9 blunt cut at PAM_start - 3 (0-based coordinate)
        cut_pos = pam_start - 3
        candidates.append(PAMCandidate(
            position=proto_start,
            strand="+",
            protospacer=protospacer,
            pam=m.group(),
            cut_position=cut_pos,
            gc_percent=calculate_gc(protospacer),
            polyt_flag=_check_polyt(protospacer),
        ))

    # ---- Reverse strand ----
    # Search the forward strand for the reverse-complement of the PAM.
    # For NGG on the reverse strand, the forward strand shows CCN.
    rev_pam = reverse_complement(pam)
    rev_pam_regex = _build_pam_regex(rev_pam)
    for m in re.finditer(rev_pam_regex, seq):
        pam_start = m.start()          # position of CCN on forward strand
        proto_end = pam_start + 3 + half  # protospacer_rc starts after CCN
        if proto_end > len(seq):
            continue
        protospacer_rc = seq[pam_start + 3:proto_end]
        protospacer = reverse_complement(protospacer_rc)
        # On the reverse strand the PAM (NGG) reads CCN on the forward strand
        # at [pam_start, pam_start+3). The protospacer lies immediately to its
        # right, and Cas9 cuts 3 bp from the PAM *inside* the protospacer, i.e.
        # at forward boundary pam_start + 3 + 3 (mirror image of the forward
        # strand rule PAM_start - 3).
        cut_pos = pam_start + 3 + 3
        # Position on forward strand (0-based start of protospacer_rc)
        fwd_pos = pam_start + 3
        candidates.append(PAMCandidate(
            position=fwd_pos,
            strand="-",
            protospacer=protospacer,
            pam=m.group(),
            cut_position=cut_pos,
            gc_percent=calculate_gc(protospacer),
            polyt_flag=_check_polyt(protospacer),
        ))

    return candidates


def check_specificity(protospacer: str,
                      genome_seq: str,
                      max_mismatch: int = 3,
                      blast_available: bool = False) -> int:
    """Check specificity of *protospacer* against *genome_seq*.

    Uses ``blastn`` if available; otherwise falls back to k-mer + Hamming
    distance matching.

    Returns
    -------
    int
        Number of genomic locations with ≤ *max_mismatch* matches.
    """
    if blast_available and shutil.which("blastn"):
        return _check_specificity_blast(protospacer, genome_seq, max_mismatch)
    return _check_specificity_fallback(protospacer, genome_seq, max_mismatch)


def _check_specificity_blast(protospacer: str, genome_seq: str,
                             max_mismatch: int) -> int:
    """Use blastn subprocess for specificity check."""
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta",
                                         delete=False) as qfh, \
             tempfile.NamedTemporaryFile(mode="w", suffix=".fasta",
                                         delete=False) as dfh:
            qfh.write(f">probe\n{protospacer}\n")
            qpath = qfh.name
            dfh.write(f">genome\n{genome_seq}\n")
            dpath = dfh.name

        db_dir = os.path.dirname(dpath)
        db_name = os.path.splitext(os.path.basename(dpath))[0]
        full_db = os.path.join(db_dir, db_name)

        subprocess.run(
            ["makeblastdb", "-in", dpath, "-dbtype", "nucl", "-out", full_db],
            check=True, capture_output=True, timeout=60
        )
        result = subprocess.run(
            ["blastn", "-query", qpath, "-db", full_db,
             "-word_size", "7", "-evalue", "1000",
             "-outfmt", "6", "-max_target_seqs", str(10000)],
            check=True, capture_output=True, text=True, timeout=60
        )

        hits = 0
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            # Column 4 (0-indexed 3) is % identity; mismatches in column 5
            # But with -outfmt 6, mismatches is not always present.
            # Use alignment length / mismatches heuristic:
            alignment_len = int(fields[3])
            mismatches = alignment_len - int(fields[2]) if len(fields) > 4 else 0
            # Actually fields[2] is %identity (float), fields[3] is alignment length
            # mismatches = alignment_len - (alignment_len * pct_id / 100)  ... approximate
            # Better: compute mismatches from %identity
            if len(fields) >= 5:
                n_mid = int(round(alignment_len * float(fields[2]) / 100))
                mismatches = alignment_len - n_mid
            if mismatches <= max_mismatch:
                hits += 1

        _cleanup_blast_db(full_db, dpath, qpath)
        return hits

    except Exception:
        return _check_specificity_fallback(protospacer, genome_seq, max_mismatch)


def _cleanup_blast_db(full_db: str, db_path: str, query_path: str):
    """Remove temporary BLAST database and query files."""
    for ext in [".nhr", ".nin", ".nsd", ".nsi", ".nsq", ".nog", ".nsd"]:
        p = full_db + ext
        if os.path.exists(p):
            os.unlink(p)
    for p in [db_path, query_path]:
        if os.path.exists(p):
            os.unlink(p)


def _check_specificity_fallback(protospacer: str, genome_seq: str,
                                max_mismatch: int) -> int:
    """Fallback specificity check using Hamming distance over sliding window."""
    protospacer = protospacer.upper()
    rev_protospacer = reverse_complement(protospacer)
    genome_seq = genome_seq.upper().replace(" ", "").replace("\n", "")
    n = len(protospacer)
    hits = _hamming_search(genome_seq, protospacer, max_mismatch)
    hits += _hamming_search(genome_seq, rev_protospacer, max_mismatch)
    return hits


def _hamming_search(haystack: str, needle: str, max_mismatch: int) -> int:
    """Count windows in *haystack* within *max_mismatch* of *needle*.

    Strategy: split the needle into ``max_mismatch + 1`` contiguous blocks.
    By the pigeonhole principle, any window with at most ``max_mismatch``
    mismatches must match at least one block exactly. We locate each block
    with ``str.find`` (C speed), collect candidate windows, then verify each
    with an early-exit Hamming check.

    This is exact — no window within tolerance is missed.
    """
    nlen = len(needle)
    hlen = len(haystack)
    if nlen == 0 or hlen < nlen:
        return 0

    n_blocks = max_mismatch + 1
    if n_blocks > nlen:
        n_blocks = nlen
    block_len = max(1, nlen // n_blocks)

    count = 0
    seen: set = set()

    for b in range(n_blocks):
        bstart = b * block_len
        bend = nlen if b == n_blocks - 1 else bstart + block_len
        block = needle[bstart:bend]
        if not block:
            continue
        start = 0
        while True:
            idx = haystack.find(block, start)
            if idx == -1:
                break
            win_start = idx - bstart
            # Window must be fully inside haystack
            if 0 <= win_start and win_start + nlen <= hlen and win_start not in seen:
                seen.add(win_start)
                mism = 0
                for a, c in zip(haystack[win_start:win_start + nlen], needle):
                    if a != c:
                        mism += 1
                        if mism > max_mismatch:
                            break
                if mism <= max_mismatch:
                    count += 1
            start = idx + 1

    return count



def rank_sgRNAs(candidates: list[PAMCandidate],
                genome_seq: str,
                config: PipelineConfig) -> list[PAMCandidate]:
    """Rank sgRNA candidates by specificity and structural quality.

    Scoring (higher is better):
    - Base score 100
    - -15 penalty for poly-T (Pol III termination signal) — in-vivo mode only;
      in-vitro guides are T7-transcribed, where a TTTT run is harmless
    - -5 penalty per off-target beyond the perfect match
    - +5 bonus if GC% in 40–65%, -10 if GC% > 75
    """
    scored: list[tuple[PAMCandidate, float]] = []
    for c in candidates:
        score = 100.0
        if not c.protospacer:
            scored.append((c, -1))
            continue
        off_targets = check_specificity(
            c.protospacer, genome_seq, config.max_mismatches,
            config.blast_available or False,
        )
        score -= 5 * max(0, off_targets - 1)
        if c.polyt_flag and config.mode == "in-vivo":
            score -= 15
        if config.mode == "in-vitro" and not c.protospacer.startswith("G"):
            c.warnings.append(
                "no 5' G: prepend a G to the protospacer for T7 transcription")
        if 40 <= c.gc_percent <= 65:
            score += 5
        elif c.gc_percent > 75:
            score -= 10
        scored.append((c, score))

    scored.sort(key=lambda x: x[1], reverse=True)
    return [c for c, _ in scored]


def design_sgRNAs(cluster_seq: str,
                  upstream_seq: str,
                  downstream_seq: str,
                  genome_seq: Optional[str] = None,
                  config: Optional[PipelineConfig] = None,
                  upstream_start: int = 0,
                  rejected_out: Optional[dict] = None) -> dict:
    """Design ranked sgRNAs for both ends of a BGC.

    Parameters
    ----------
    cluster_seq : str
        BGC sequence.
    upstream_seq : str
        Sequence upstream of the BGC (5' flank).
    downstream_seq : str
        Sequence downstream of the BGC (3' flank).
    genome_seq : str, optional
        Full genome for specificity filtering. If ``None``, the concatenation of
        upstream + cluster + downstream is used.
    config : PipelineConfig, optional
        Configuration overrides (defaults if omitted).
    upstream_start : int
        Coordinate, in ``genome_seq``, of the first base of ``upstream_seq``
        (``cluster.start - len(upstream_seq)``). Cut positions are reported in
        ``genome_seq`` coordinates, so this MUST be set whenever ``genome_seq``
        is a larger sequence than upstream + cluster + downstream.
    rejected_out : dict, optional
        If given, filled with ``{"left": [...], "right": [...]}`` listing
        candidates dropped for cleaving inside the BGC (with their sites).

    Returns
    -------
    dict
        ``{"left": [...], "right": [...]}`` — top-N ranked candidates per end.
    """
    if config is None:
        config = PipelineConfig()
    if genome_seq is None:
        genome_seq = upstream_seq + cluster_seq + downstream_seq

    left_candidates = find_pam_sites(upstream_seq, config.pam_sequence,
                                     config.protospacer_length)
    right_candidates = find_pam_sites(downstream_seq, config.pam_sequence,
                                      config.protospacer_length)

    # Offset cut positions to be relative to genome_seq
    upstream_offset = upstream_start
    cluster_start = upstream_start + len(upstream_seq)
    cluster_end = cluster_start + len(cluster_seq)
    downstream_offset = cluster_end
    left_candidates = _offset_candidates(left_candidates, upstream_offset)
    right_candidates = _offset_candidates(right_candidates, downstream_offset)

    left_ranked = rank_sgRNAs(left_candidates, genome_seq, config)
    right_ranked = rank_sgRNAs(right_candidates, genome_seq, config)

    # Drop guides that would also cleave inside the BGC. One index over the
    # local window (upstream + cluster + downstream) serves every candidate;
    # pair-level checks over the full fragment happen later.
    local = upstream_seq + cluster_seq + downstream_seq
    idx = build_site_index(local, offset=upstream_start,
                           protospacer_len=config.protospacer_length,
                           seed_len=config.seed_length)
    kept = {}
    rejected = {"left": [], "right": []}
    for side, ranked in (("left", left_ranked), ("right", right_ranked)):
        kept[side] = []
        for c in ranked:
            sites = find_cut_sites(c.protospacer, idx, config.max_mismatches)
            if side == "left":     # fragment starts at the cut: (cut, cluster_end]
                inside = [x for x in sites if c.cut_position < x.cut_position <= cluster_end]
            else:                  # fragment ends at the cut: [cluster_start, cut)
                inside = [x for x in sites if cluster_start <= x.cut_position < c.cut_position]
            high = [x for x in inside if x.severity == "high"]
            if high and config.exclude_internal_cuts:
                rejected[side].append((c, inside))
                continue
            c.internal_cuts = inside
            for x in inside:
                c.warnings.append(
                    f"possible additional cut inside the BGC at {x.cut_position} "
                    f"(PAM {x.pam}, {x.mismatches} mismatch(es), severity {x.severity})")
            kept[side].append(c)
    if rejected_out is not None:
        rejected_out.update(rejected)

    return {
        "left": kept["left"][:config.top_n_sgRNAs],
        "right": kept["right"][:config.top_n_sgRNAs],
    }


def _offset_candidates(candidates: list[PAMCandidate],
                       offset: int) -> list[PAMCandidate]:
    """Adjust position and cut_position of each candidate by offset."""
    adjusted = []
    for c in candidates:
        adjusted.append(dataclasses.replace(
            c, position=c.position + offset, cut_position=c.cut_position + offset))
    return adjusted


__all__ = [
    "find_pam_sites",
    "calculate_gc",
    "check_specificity",
    "rank_sgRNAs",
    "design_sgRNAs",
    "complement",
    "reverse_complement",
    "_offset_candidates",
]
