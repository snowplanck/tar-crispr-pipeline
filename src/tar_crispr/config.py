from dataclasses import dataclass, field
from typing import Optional

@dataclass
class PipelineConfig:
    pam_window: int = 500
    protospacer_length: int = 20
    pam_sequence: str = "NGG"
    top_n_sgRNAs: int = 5
    homology_arm_length: int = 50
    homology_arm_range: tuple = (40, 60)
    max_shift: int = 100
    shift_increment: int = 5
    min_gc: float = 0.40
    max_gc: float = 0.65
    max_gc_warning: float = 0.75
    polyt_threshold: int = 4
    primer_tm_min: float = 58.0
    primer_tm_max: float = 62.0
    primer_max_delta_tm: float = 2.0
    primer_anneal_min: int = 18
    primer_anneal_max: int = 25
    max_mismatches: int = 3
    avoid_enzymes: list = field(default_factory=list)
    # Cas9 delivery mode. "in-vitro": HMW genomic DNA is digested with
    # Cas9/sgRNA RNP before transformation (vector never sees Cas9).
    # "in-vivo": Cas9 is expressed in yeast, so the vector and the yeast
    # genome are also potential cleavage targets.
    mode: str = "in-vitro"
    seed_length: int = 12
    exclude_internal_cuts: bool = True
    # "mit": genome-wide, position-weighted off-target score (Hsu 2013-style);
    # "legacy": plain count of look-alike sites (old behaviour).
    specificity_model: str = "mit"
    n_pairs: int = 5          # guide pairs for which homology arms are evaluated
    min_specificity: float = 50.0   # warn below this guide specificity (0-100)
    blast_available: Optional[bool] = None
    rnafold_available: Optional[bool] = None

@dataclass
class ClusterInfo:
    start: int
    end: int
    name: str = "BGC"
    description: str = ""
    strand: str = "+"
    gc_percent: Optional[float] = None

    @property
    def length(self) -> int:
        """Length of the cluster in bp (end - start)."""
        return self.end - self.start

@dataclass
class PAMCandidate:
    position: int
    strand: str
    protospacer: str
    pam: str
    cut_position: int
    gc_percent: float
    polyt_flag: bool
    internal_cuts: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    specificity_score: Optional[float] = None   # 0-100, higher = more specific
    n_offtargets: int = 0

@dataclass
class CutSite:
    """A position where a guide can plausibly direct Cas9 cleavage."""
    cut_position: int          # 0-based boundary (same convention as PAMCandidate)
    strand: str
    pam: str
    pam_class: str             # "NGG", "NAG" or "NGA"
    mismatches: int            # total protospacer mismatches
    seed_mismatches: int       # mismatches in the PAM-proximal seed
    severity: str              # "high" or "medium"
    locus: str = "fragment"    # "fragment", "vector", "yeast" or "genome"
    site: str = ""             # genomic protospacer sequence at the site
    score: float = 0.0         # off-target hit score (0-1), genome scan only
    mismatch_positions: tuple = ()  # 0 = PAM-distal ... 19 = PAM-proximal

@dataclass
class HomologyArm:
    sequence: str
    end: str
    length: int
    gc_percent: float
    start_pos: int
    end_pos: int
    uniqueness: bool
    secondary_structure: bool
    issues: list = field(default_factory=list)

@dataclass
class TailedPrimer:
    name: str
    sequence: str
    tail: str
    annealing_region: str
    tm: float
    gc_percent: float
    hairpin: bool
    self_dimer: bool = False
    cross_dimer: bool = False
    issues: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    valid: bool = True

@dataclass
class AssemblyResult:
    success: bool
    final_size: int
    circular: bool
    issues: list = field(default_factory=list)
    final_sequence: str = ""

__all__ = [
    "PipelineConfig",
    "ClusterInfo",
    "PAMCandidate",
    "CutSite",
    "HomologyArm",
    "TailedPrimer",
    "AssemblyResult",
]
