"""Step 3: Definition of excised fragment ends after Cas9 cleavage.

SpCas9 produces a blunt double-strand break 3 bp upstream of the PAM.
Given one sgRNA per end, this module computes the exact linear fragment
that remains — from the left cut junction to the right cut junction —
containing the full BGC.
"""
from dataclasses import dataclass
from typing import Optional
from Bio.SeqUtils import gc_fraction
from tar_crispr.config import PAMCandidate


@dataclass
class FragmentEnds:
    """Represents the fragment resulting from Cas9 cutting at both ends."""
    sequence: str
    left_sgRNA: PAMCandidate
    right_sgRNA: PAMCandidate
    left_cut_pos: int
    right_cut_pos: int
    fragment_length: int
    gc_percent: float
    left_cut_sequence: str
    right_cut_sequence: str


def compute_cut_position(sgRNA: PAMCandidate, strand_context: str) -> int:
    """Compute the exact 0-based cut position for a given sgRNA.

    Parameters
    ----------
    sgRNA : PAMCandidate
        The chosen sgRNA candidate.
    strand_context : str
        ``'left'`` or ``'right'`` — which end of the BGC this sgRNA targets.

    Returns
    -------
    int
        0-based coordinate of the cut site in the full genomic sequence.
    """
    return sgRNA.cut_position


def extract_fragment(full_sequence: str,
                     left_sgRNA: PAMCandidate,
                     right_sgRNA: PAMCandidate) -> FragmentEnds:
    """Extract the linear fragment between two Cas9 cut sites.

    Parameters
    ----------
    full_sequence : str
        The full genomic sequence (upstream flank + BGC + downstream flank).
    left_sgRNA : PAMCandidate
        sgRNA for the left (5') cut.
    right_sgRNA : PAMCandidate
        sgRNA for the right (3') cut.

    Returns
    -------
    FragmentEnds
        Dataclass with the fragment sequence and metadata.
    """
    left_cut = left_sgRNA.cut_position
    right_cut = right_sgRNA.cut_position

    if left_cut >= right_cut:
        raise ValueError(
            f"Left cut ({left_cut}) must be before right cut ({right_cut}). "
            "Check that upstream sgRNA precedes downstream sgRNA."
        )

    fragment_seq = full_sequence[left_cut:right_cut]

    if len(fragment_seq) == 0:
        raise ValueError("Fragment is empty — cut positions too close.")

    gc = round(gc_fraction(fragment_seq) * 100, 2) if len(fragment_seq) > 0 else 0.0

    return FragmentEnds(
        sequence=fragment_seq,
        left_sgRNA=left_sgRNA,
        right_sgRNA=right_sgRNA,
        left_cut_pos=left_cut,
        right_cut_pos=right_cut,
        fragment_length=len(fragment_seq),
        gc_percent=gc,
        left_cut_sequence=fragment_seq[:20],
        right_cut_sequence=fragment_seq[-20:],
    )


def get_cut_junction_sequences(sgRNA: PAMCandidate, full_sequence: str,
                               num_nt: int = 20) -> tuple[str, str]:
    """Extract sequences flanking the cut site on both strands.

    Parameters
    ----------
    sgRNA : PAMCandidate
        The sgRNA whose cut site we examine.
    full_sequence : str
        Full genomic sequence.
    num_nt : int
        Number of nucleotides to extract on each side of the cut.

    Returns
    -------
    tuple of (upstream_flank, downstream_flank)
        Sequences immediately 5' and 3' of the cut site (on the forward strand).
    """
    cut = sgRNA.cut_position
    upstream = full_sequence[max(0, cut - num_nt):cut]
    downstream = full_sequence[cut:cut + num_nt]
    return upstream, downstream


__all__ = [
    "FragmentEnds",
    "compute_cut_position",
    "extract_fragment",
    "get_cut_junction_sequences",
]
