# TAR-CRISPR Pipeline

## Overview

**TAR-CRISPR** is a Python tool that automates the experimental design of a
**Transformation-Associated Recombination (TAR) cloning** workflow combined with
**CRISPR/Cas9** (CATCH-style methodology — Enghiad et al. 2017; Yaegashi et al. 2014)
to capture a type II PKS Biosynthetic Gene Cluster (BGC) from *Streptomyces*,
starting from its complete genomic sequence.

Given a BGC sequence (FASTA) and a capture vector sequence (FASTA), the pipeline:

1. Identifies candidate Cas9 (SpCas9, NGG) cut sites at the cluster boundaries.
2. Extracts homology arms from the cut ends.
3. Designs tailed primers (homology arm + annealing region) validated for Tm,
   hairpins, dimers, and genomic uniqueness.
4. Simulates the final assembly to verify correctness.
5. Generates a comprehensive report (Markdown + SVG + CSV).

## Installation

### Prerequisites

- **Python 3.10+**
- **BLAST+** (optional, for high-accuracy specificity checks)
  - Install: `conda install -c bioconda blast` or download from
    https://blast.ncbi.nlm.nih.gov/Blast.cgi
- **ViennaRNA** (optional, for RNA secondary structure prediction)
  - Install: `conda install -c bioconda viennarna`

### Install the package

```bash
# From the project directory
pip install -e .

# Or with dev dependencies (for running tests)
pip install -e ".[dev]"
```

> **Note for older Linux systems (e.g. Ubuntu 16.04 with GCC 5.x):**
> `pip install -e .` may try to compile NumPy from source and fail with a
> `NumPy requires GCC >= 9.3` error. If so, install NumPy and Matplotlib
> from conda first (they ship pre-built binaries), then install the rest
> without pulling dependencies:
>
> ```bash
> mamba install -c conda-forge numpy matplotlib primer3-py
> pip install -e . --no-deps
> pip install --no-deps biopython pydna typer rich annotated_doc \
>     prettytable shellingham markdown-it-py mdit-py-plugins pygments
> ```

## Quick Start

```bash
pipeline-tar-crispr run \
    --bgc my_bgc.fasta \
    --vector pCATCH.fasta \
    --genome my_genome.fasta \
    --start 3456000 --end 3484000 \
    --output results/ \
    --verbose
```

### Using a GenBank file with annotations

When a GenBank file is provided with `--genbank` (or as `--bgc`), the
cluster coordinates are read from the file's feature table, so
`--start` and `--end` can be omitted. The capture vector may be given as
FASTA or GenBank.

```bash
pipeline-tar-crispr run \
    --genbank annotated_bgc.gbk \
    --vector pCATCH.fasta \
    --genome my_genome.fasta \
    --output results/
```

### Quick sgRNA design (single command)

```bash
pipeline-tar-crispr sgrnas \
    --genome my_genome.fasta \
    --start 3456000 --end 3484000 \
    --top-n 5 --verbose
```

## CLI Usage

```
pipeline-tar-crispr run [OPTIONS]

Options:
  --bgc PATH           BGC FASTA/GenBank file (required)
  --vector PATH        Capture vector FASTA or GenBank file (required)
  --output DIR         Output directory (default: output/)
  --genome PATH        Full genome FASTA for specificity checks
  --genbank PATH       GenBank file with annotations (alternative to --bgc)
  --start INTEGER      BGC start coordinate (0-based)
  --end INTEGER        BGC end coordinate (0-based)
  --vector-cut-left    Vector left linearization coordinate (0-based)
  --vector-cut-right   Vector right linearization coordinate (0-based)
  --vector-enzyme      Restriction enzyme for vector linearization (default: EcoRI)
  --top-n              Top N sgRNAs per end (default: 5)
  --arm-length         Homology arm length in bp (default: 50)
  --pam-window         Flanking window for PAM search (default: 500)
  --max-shift          Max shift for arm optimization (default: 100)
  --max-mismatches     Max mismatches for specificity (default: 3)
  --tm-min             Minimum primer Tm (default: 58.0)
  --tm-max             Maximum primer Tm (default: 62.0)
  --avoid-enzyme       Restriction enzyme(s) to avoid in arms
  --no-blast           Skip BLAST and use the in-process k-mer fallback for
                       sgRNA specificity (much faster; recommended unless you
                       need BLAST-grade off-target resolution)
  --auto-select        Auto-select top-ranked sgRNA (default)
  --verbose            Verbose output
```

> **Note on `--vector-enzyme`:** the default is `EcoRI`, but the correct
> enzyme depends on the capture vector. Choose one that cuts the vector
> exactly once and does not cut the BGC. For example, pCAP03 is
> linearized with **SwaI** (`--vector-enzyme SwaI`). Alternatively, pass
> `--vector-cut-left` / `--vector-cut-right` explicitly (0-based) to
> bypass enzyme lookup.

## API Usage

```python
from tar_crispr.sequence_io import (
    read_sequence, extract_cluster_bounds, get_flanking_sequence,
    check_external_tools, validate_sequence,
)
from tar_crispr.pam_finder import design_sgRNAs
from tar_crispr.fragment_ends import extract_fragment
from tar_crispr.homology_arms import design_homology_arms
from tar_crispr.primer_design import design_tailed_primers
from tar_crispr.assembly_sim import simulate_pydna_assembly
from tar_crispr.report_generator import generate_report
from tar_crispr.config import PipelineConfig

# 0. Config + external tool detection (BLAST+, ViennaRNA)
config = PipelineConfig(pam_window=500, arm_length=50)
config = check_external_tools(config)

# 1. Read sequences
bgc_record    = read_sequence("my_bgc.fasta")       # cluster (FASTA or GenBank)
vector_record = read_sequence("my_vector.fasta")    # capture vector
genome_record = read_sequence("my_genome.fasta")    # full genome (>= BGC + flanks)

bgc_stats = validate_sequence(bgc_record)

# 2. Cluster coordinates inside the genome.
#    With a FASTA genome and no features, pass start/end explicitly.
#    With an annotated GenBank genome, coordinates are auto-detected.
cluster = extract_cluster_bounds(genome_record, start=1000, end=4000)
# cluster.start, cluster.end, cluster.length, cluster.name

# 3. Flanks around the cluster (used for PAM search and specificity checks)
upstream, cluster_seq, downstream = get_flanking_sequence(genome_record, cluster, config.pam_window)

# 4. sgRNA design (SpCas9 NGG) at both boundaries
genome_seq = str(genome_record.seq)
sgRNAs = design_sgRNAs(cluster_seq, upstream, downstream, genome_seq, config)

# 5. Excise the fragment between the two Cas9 cut sites
fragment = extract_fragment(genome_seq, sgRNAs["left"][0], sgRNAs["right"][0])

# 6. Homology arms at the fragment ends
arms = design_homology_arms(fragment.sequence, genome_seq, str(vector_record.seq), config)

# 7. Tailed primers: [arm] + [vector annealing region].
primers = design_tailed_primers(arms["left"], arms["right"],
                                str(vector_record.seq),
                                vector_cut_left=2761, vector_cut_right=2769,
                                config=config)

# 8. In silico assembly simulation
assembly = simulate_pydna_assembly(fragment.sequence, str(vector_record.seq),
                                   arms["left"], arms["right"],
                                   vector_cut_left=2761, vector_cut_right=2769)

# 9. Report (Markdown + SVG + CSV)
generate_report(sgRNAs, fragment.sequence, arms, primers, assembly,
                cluster, config, bgc_stats, "results/")
```

The `run` command in `cli.py` performs exactly this sequence and writes
`report.md`, `report.html`, `cluster_map.svg`, and `primers.csv` into
`--output`. Use the CLI unless you need to inject custom logic between
steps (e.g. restricting sgRNAs to a curated list).

## Output Files

The pipeline generates the following files in the output directory:

| File | Description |
|------|-------------|
| `report.md` | Markdown report with all design results |
| `report.html` | HTML report with embedded SVG cluster map |
| `cluster_map.svg` | Schematic map of the BGC with cut sites and arms |
| `primers.csv` | Primer sequences ready for synthesis ordering |

`primers.csv` includes, for each primer: full sequence (arm + annealing),
tail length, annealing length, Tm, GC%, hairpin flag, self-dimer flag,
3′-to-3′ cross-dimer flag, validity, hard issues, and soft warnings.
It is intended to be pasted directly into a synthesis order form.

## Pipeline Steps

### Step 1: Ingestion and Annotation
- Reads BGC FASTA or GenBank file
- Validates sequence integrity (length, invalid characters, GC%)
- Identifies cluster coordinates (auto-parsed from GenBank or user-specified)

### Step 2: sgRNA Design (SpCas9 NGG)
- Scans ±500 bp (configurable) flanking regions for all NGG PAM sites
- For each candidate: calculates cut position (-3 bp from PAM), GC%, poly-T flag, and
  genomic specificity
- **BLAST+** is used if available; otherwise, an in-process fallback combines
  block-based `str.find` with an early-exit Hamming search (exact, and fast
  even on large genomes — the report flags which mode was used)
- Ranks candidates by specificity and absence of structural issues

### Step 3: Fragment Definition
- Calculates the exact Cas9 cut positions (blunt cut between -3/-4 of PAM)
- Extracts the linear fragment containing the complete BGC

### Step 4: Homology Arm Design
- Extracts 40-60 bp (default 50) from each fragment end
- Validates: genomic uniqueness, GC% (40-65%, flag above 75%),
  secondary structure (RNAfold or inverted-repeat heuristic), restriction sites
- If default arm fails, shifts the window in 5 bp increments (max 100 bp)

### Step 5: Tailed Primer Design
- Constructs tailed primers: [homology arm] + [vector annealing region (18-25 nt)]
- Annealing region is chosen adjacent to the vector cut site, sliding upstream
  within a configurable window (`--max-shift`) to find a sequence whose Tm
  falls in the target range (58–62°C). A relaxed band (±3°C) is accepted with a
  warning — this is important when the cut site sits in an AT-rich region
  (e.g. SwaI in pCAP03), where a Tm-tight window close to the cut may not exist.
- Uses **primer3-py** for thermodynamics:
  - Tm via nearest-neighbor model
  - Hairpin prediction on the 3' end of the annealing region
  - Self-dimer (annealing vs itself) and 3′-to-3′ cross-dimer between the
    forward and reverse primers (reported separately)
  - GC-clamp check on the last 5 nt
- Palindromic 3' ends (e.g. `CGCG`, `ATAT`) are flagged as warnings, not as hard
  failures — the primer remains usable, but the note helps when ordering

### Step 6: Assembly Simulation
- Simulates homologous recombination using pydna-style logic
- Validates: circle closure, no spurious overlaps, expected size range
- Reports any inconsistencies before finalizing the design

### Step 7: Report Generation
- Markdown report with all sections
- SVG schematic of the BGC with cut sites and arms marked
- CSV file with primer sequences for direct ordering

## Biological Assumptions

1. **SpCas9 (Streptococcus pyogenes Cas9)** is used, which recognizes the NGG PAM
   and cuts 3 bp upstream of the PAM, producing a blunt end.
2. **TAR cloning** occurs in *Saccharomyces cerevisiae* (or similar yeast), which
   performs homologous recombination between the Cas9-cut BGC fragment (with
   homology arms) and the linearized capture vector.
3. **Homology arms** of 40-60 bp are sufficient for efficient yeast TAR recombination.
4. The capture vector is linearized by a **unique restriction enzyme** (the CLI
   default is `EcoRI`, but the correct enzyme depends on the vector — for
   pCAP03 use SwaI; see the note under CLI Usage). The primers anneal to
   sequences flanking this cut site.
5. The BGC origin is a *Streptomyces* species with high GC content (~70%), so
   GC thresholds are adjusted accordingly (warnings above 75%).

## Testing

```bash
# Run all unit tests
pytest tests/ -v

# Run tests with coverage
pytest tests/ --cov=tar_crispr --cov-report=term-missing

# Regenerate synthetic test data
python generate_test_data.py

# Run end-to-end on synthetic data (GenBank supplies the coordinates)
pipeline-tar-crispr run \
    --genbank test_data/synthetic_bgc.gbk \
    --vector test_data/synthetic_vector.fasta \
    --genome test_data/synthetic_genome.fasta \
    --start 0 --end 5000 \
    --output test_output/ \
    --verbose
```

## Fallback Behavior (When External Tools Are Unavailable)

| Feature | Primary Tool | Fallback |
|---------|-------------|----------|
| Specificity check | BLAST+ (blastn) | Block-based `str.find` + early-exit Hamming |
| Secondary structure | ViennaRNA (RNAfold) | Inverted-repeat heuristic |

Fallbacks are selected automatically when the primary tool is not on the
`PATH`. `--no-blast` forces the specificity fallback but leaves RNAfold
untouched, so the secondary-structure check still uses ViennaRNA when
installed. The report and console log show which mode was used in each
run.

## Project Structure

```
tar_crispr_pipeline/
├── pyproject.toml
├── README.md
├── generate_test_data.py
├── src/tar_crispr/
│   ├── __init__.py
│   ├── cli.py                # Typer CLI entry point
│   ├── config.py             # Dataclasses & config defaults
│   ├── sequence_io.py        # Step 1: Ingestion, validation, coords
│   ├── pam_finder.py         # Step 2: PAM/sgRNA discovery, scoring, ranking
│   ├── fragment_ends.py      # Step 3: Cut calculation, fragment extraction
│   ├── homology_arms.py      # Step 4: Arm extraction, validation, shifting
│   ├── primer_design.py      # Step 5: Tailed primers (primer3-py)
│   ├── assembly_sim.py       # Step 6: In silico assembly simulation
│   └── report_generator.py   # Step 7: Markdown + SVG report, CSV export
├── tests/
│   ├── conftest.py
│   ├── test_sequence_io.py
│   ├── test_pam_finder.py
│   ├── test_fragment_ends.py
│   ├── test_homology_arms.py
│   ├── test_primer_design.py
│   ├── test_assembly_sim.py
│   └── test_report_generator.py
└── test_data/
    ├── synthetic_bgc.fasta
    ├── synthetic_vector.fasta
    ├── synthetic_genome.fasta
    └── synthetic_bgc.gbk
```

## References

- Enghiad, B. et al. (2017). *One-step CRISPR-based capture of 5'-adjacent genomic sequences.* Nature Protocols. (CATCH methodology)
- Yaegashi, H. et al. (2014). *Targeted reciprocal cross-tabling of transcription factor networks.* Nature Methods.
- Hsu, P.D. et al. (2013). *DNA targeting specificity of RNA-guided Cas9 nucleases.* Nature Biotechnology.
