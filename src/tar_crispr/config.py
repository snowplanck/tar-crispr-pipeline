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
    blast_available: Optional[bool] = None
    rnafold_available: Optional[bool] = None

@dataclass
class ClusterInfo:
    start: int
    end: int
    name: str = "BGC"
    description: str = ""
    strand: str = "+"

@dataclass
class PAMCandidate:
    position: int
    strand: str
    protospacer: str
    pam: str
    cut_position: int
    gc_percent: float
    polyt_flag: bool

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
    dimer: bool
    issues: list = field(default_factory=list)
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
    "HomologyArm",
    "TailedPrimer",
    "AssemblyResult",
]
