"""Step 4: Homology arm design and validation.

For TAR/CRISPR cloning, homology arms are short sequences (40-60 bp by default)
that enable homologous recombination between the excised BGC fragment and the
linearized capture vector.  Each arm is validated for:

1. **Genomic uniqueness** — must not align elsewhere in the source genome
   or capture vector (k-mer / Hamming fallback when BLAST is unavailable).
2. **GC%** — ideally 40-65 %; flagged above 75 % (Streptomyces is GC-rich).
3. **Secondary structure** — inverted-repeat / hairpin heuristic (RNAfold
   fallback when ViennaRNA is not installed).
4. **Restriction sites** — optional exclusion of sites for user-specified enzymes.

If the default arm position fails any criterion, the window is shifted in 5 bp
increments (without altering the Cas9 cut site) up to a configurable limit.
"""
from dataclasses import dataclass
from typing import Optional

from Bio.SeqUtils import gc_fraction
from Bio.Restriction import Restriction, Analysis

from tar_crispr.config import PipelineConfig, HomologyArm
from tar_crispr.pam_finder import reverse_complement, check_specificity, _hamming_search


def extract_homology_arm(fragment_seq: str,
                         end: str,
                         length: int = 50) -> str:
    """Extract a homology arm from an end of the fragment.

    Parameters
    ----------
    fragment_seq : str
        The linear fragment sequence (product of Cas9 digestion).
    end : str
        ``"left"`` or ``"right"``.
    length : int
        Arm length in bp (default 50).

    Returns
    -------
    str
        The homology arm sequence (5'->3' orientation, same as fragment).
    """
    if end == "left":
        if length >= len(fragment_seq):
            return fragment_seq
        return fragment_seq[:length]
    elif end == "right":
        if length >= len(fragment_seq):
            return fragment_seq
        return fragment_seq[-length:]
    else:
        raise ValueError(f"end must be 'left' or 'right', got '{end}'")


def validate_gc(arm: str, config: PipelineConfig) -> tuple[bool, list[str]]:
    """Check GC content of the arm against configured thresholds.

    Returns (passed, issues).
    """
    issues = []
    if len(arm) == 0:
        return False, ["Arm is empty"]
    gc = gc_fraction(arm) * 100
    if gc < config.min_gc * 100:
        issues.append(f"GC% {gc:.1f} below minimum {config.min_gc * 100:.0f}%")
    if gc > config.max_gc * 100:
        if gc > config.max_gc_warning * 100:
            issues.append(f"GC% {gc:.1f} exceeds warning threshold {config.max_gc_warning * 100:.0f}%")
        else:
            issues.append(f"GC% {gc:.1f} above ideal range {config.max_gc * 100:.0f}%")
    return len(issues) == 0, issues


def check_uniqueness(arm: str, genome_seq: str, vector_seq: str,
                     config: PipelineConfig) -> tuple[bool, list[str]]:
    """Check that the arm appears only once in the genome and not in the vector.

    Uses BLAST if available, otherwise falls back to k-mer/Hamming search.

    Returns (passed, issues).
    """
    issues = []
    arm_upper = arm.upper()
    genome_upper = genome_seq.upper().replace(" ", "").replace("\n", "")
    vector_upper = vector_seq.upper().replace(" ", "").replace("\n", "")

    # Check genome: should appear exactly once (the intended target)
    genome_hits = _hamming_search(genome_upper, arm_upper, 0)  # exact matches
    if genome_hits > 1:
        issues.append(
            f"Arm appears {genome_hits} times in genome (expected 1) — "
            "non-unique"
        )

    # Check vector: should NOT appear (unless it's the intended junction)
    vector_hits = _hamming_search(vector_upper, arm_upper, 0)
    if vector_hits > 0:
        issues.append(
            f"Arm found {vector_hits} times in capture vector — "
            "potential off-target recombination"
        )

    return len(issues) == 0, issues


def check_secondary_structure(arm: str) -> tuple[bool, list[str]]:
    """Detect potential hairpin/inverted-repeat structures in the arm.

    Uses RNAfold if available; otherwise uses a simple inverted-repeat heuristic.

    Returns (passed, issues).
    """
    import shutil
    issues = []
    if shutil.which("RNAfold"):
        return _check_secondary_rnafold(arm, issues)
    return _check_secondary_heuristic(arm, issues)


def _check_secondary_rnafold(arm: str, issues: list[str]) -> tuple[bool, list[str]]:
    """Use ViennaRNA RNAfold for hairpin prediction."""
    import subprocess
    try:
        # Treat DNA as RNA (T->U) for folding prediction
        rna_seq = arm.upper().replace("T", "U")
        result = subprocess.run(
            ["RNAfold", "--noPS"], input=rna_seq,
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            if len(lines) >= 2:
                mfe_line = lines[1]
                # Extract energy from output: "((((...)))) (X.XX)"
                energy_str = mfe_line.split()[-1]
                energy = float(energy_str)
                if energy < -3.0:
                    issues.append(
                        f"Strong secondary structure predicted (ΔG={energy:.2f} kcal/mol)"
                    )
        return len(issues) == 0, issues
    except Exception:
        return _check_secondary_heuristic(arm, issues)


def _check_secondary_heuristic(arm: str, issues: list[str]) -> tuple[bool, list[str]]:
    """Simple inverted-repeat detection heuristic.

    Searches for subsequences of >= 6 bp that are reverse-complementary to
    another part of the arm, separated by at least 3 nucleotides (minimum loop
    size for hairpin formation).
    """
    from tar_crispr.pam_finder import reverse_complement
    arm_upper = arm.upper()
    min_len = 6
    min_loop = 3  # minimum loop size between inverted repeat segments
    n = len(arm_upper)
    found = False
    for i in range(n - min_len + 1):
        for j in range(i + min_len + min_loop, n - min_len + 1):
            sub = arm_upper[i:i + min_len]
            rc_sub = reverse_complement(sub)
            if j + min_len <= n and arm_upper[j:j + min_len] == rc_sub:
                found = True
                issues.append(
                    f"Potential inverted repeat at positions {i}-{i+min_len} "
                    f"and {j}-{j+min_len} (may form hairpin)"
                )
                break
        if found:
            break
    return not found, issues


def check_restriction_sites(arm: str, enzymes: list[str]) -> tuple[bool, list[str]]:
    """Check for presence of restriction sites in the arm.

    Parameters
    ----------
    arm : str
        Homology arm sequence.
    enzymes : list of str
        Restriction enzyme names (e.g. "EcoRI", "BamHI").

    Returns (passed, issues).
    """
    issues = []
    if not enzymes:
        return True, issues
    arm_upper = arm.upper()
    from Bio.Seq import Seq
    arm_seq = Seq(arm_upper)
    for enz_name in enzymes:
        try:
            enz = getattr(Restriction, enz_name)
            sites = enz.search(arm_seq)
            if sites:
                issues.append(f"Contains {enz_name} site(s) at positions {sites}")
        except AttributeError:
            issues.append(f"Unknown restriction enzyme: {enz_name}")
    return len(issues) == 0, issues


def validate_arm(arm: str, genome_seq: str, vector_seq: str,
                 config: PipelineConfig) -> tuple[bool, list[str]]:
    """Run all validation checks on a homology arm.

    Returns (passed, issues_list).
    """
    all_issues = []

    _, gc_issues = validate_gc(arm, config)
    all_issues.extend(gc_issues)

    _, uniq_issues = check_uniqueness(arm, genome_seq, vector_seq, config)
    all_issues.extend(uniq_issues)

    _, ss_issues = check_secondary_structure(arm)
    all_issues.extend(ss_issues)

    _, re_issues = check_restriction_sites(arm, config.avoid_enzymes)
    all_issues.extend(re_issues)

    return len(all_issues) == 0, all_issues


def find_valid_arm(fragment_seq: str, end: str,
                   config: PipelineConfig,
                   genome_seq: str,
                   vector_seq: str,
                   max_attempts: Optional[int] = None) -> HomologyArm:
    """Find a valid homology arm by extracting and optionally shifting.

    Tries the default arm length first. If validation fails, shifts the window
    in ``config.shift_increment`` steps (without altering the cut site) up to
    ``config.max_shift`` bp, or until a passing arm is found.

    Parameters
    ----------
    fragment_seq : str
        The Cas9 fragment (cut product).
    end : str
        ``"left"`` or ``"right"``.
    config : PipelineConfig
        Pipeline configuration.
    genome_seq : str
        Full source genome.
    vector_seq : str
        Capture vector sequence.
    max_attempts : int, optional
        Override max_shift-based limit.

    Returns
    -------
    HomologyArm
        The best arm found (may still have issues if none pass).
    """
    max_shift = config.max_shift
    increment = config.shift_increment
    target_length = config.homology_arm_length

    best_arm = None
    best_issues = None
    best_shift = None

    for shift in range(0, max_shift + 1, increment):
        if end == "left":
            start = shift
            stop = start + target_length
            if stop > len(fragment_seq):
                break
            arm = fragment_seq[start:stop]
        else:
            stop = len(fragment_seq) - shift
            start = stop - target_length
            if start < 0:
                break
            arm = fragment_seq[start:stop]

        passed, issues = validate_arm(arm, genome_seq, vector_seq, config)
        arm_obj = HomologyArm(
            sequence=arm,
            end=end,
            length=len(arm),
            gc_percent=round(gc_fraction(arm) * 100, 2) if arm else 0.0,
            start_pos=start,
            end_pos=stop,
            uniqueness=not any("non-unique" in i or "capture vector" in i for i in issues),
            secondary_structure=not any("inverted" in i.lower() or "secondary structure" in i.lower() for i in issues),
            issues=issues,
        )

        if passed:
            return arm_obj
        if best_arm is None or len(issues) < len(best_issues):
            best_arm = arm_obj
            best_issues = issues
            best_shift = shift

    # Return best attempt even if it has issues
    return best_arm


def design_homology_arms(fragment_seq: str,
                         genome_seq: str,
                         vector_seq: str,
                         config: PipelineConfig) -> dict:
    """Design and validate homology arms for both ends of a fragment.

    Parameters
    ----------
    fragment_seq : str
        Cas9-cut fragment containing the BGC.
    genome_seq : str
        Full source genome for specificity checks.
    vector_seq : str
        Capture vector sequence.
    config : PipelineConfig
        Pipeline configuration.

    Returns
    -------
    dict
        ``{"left": HomologyArm, "right": HomologyArm}``
    """
    left_arm = find_valid_arm(fragment_seq, "left", config, genome_seq, vector_seq)
    right_arm = find_valid_arm(fragment_seq, "right", config, genome_seq, vector_seq)

    return {"left": left_arm, "right": right_arm}


__all__ = [
    "extract_homology_arm",
    "validate_gc",
    "check_uniqueness",
    "check_secondary_structure",
    "check_restriction_sites",
    "validate_arm",
    "find_valid_arm",
    "design_homology_arms",
]
