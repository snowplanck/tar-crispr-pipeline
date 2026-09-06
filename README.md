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

## Quick Start

```bash
pipeline-tar-crispr run \
    --bgc my_bgc.fasta \
    --vector pCATCH.fasta \
    --genome my_genome.fasta \
    --output results/ \
    --verbose
```

### Using a GenBank file with annotations

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
  --vector PATH        Capture vector FASTA file (required)
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
  --no-blast           Force k-mer fallback (skip BLAST)
  --auto-select        Auto-select top-ranked sgRNA (default)
  --verbose            Verbose output
```

## API Usage

```python
from tar_crispr.sequence_io import read_sequence, extract_cluster_bounds, get_flanking_sequence, check_external_tools
from tar_crispr.pam_finder import design_sgRNAs
from tar_crispr.fragment_ends import extract_fragment
from tar_crispr.homology_arms import design_homology_arms
from tar_crispr.primer_design import design_tailed_primers
from tar_crispr.assembly_sim import simulate_pydna_assembly
from tar_crispr.report_generator import generate_report
from tar_crispr.config import PipelineConfig

config = PipelineConfig()
config = check_external_tools(config)

# Step 1: Read sequences
bgc_record = read_sequence("my_bgc.fasta")
vector_record = read_sequence("vector.fasta")
genome_record = read_sequence("my_genome.fasta")
cluster = extract_cluster_bounds(bgc_record, start=1000, end=4000)
upstream, cluster_seq, downstream = get_flanking_sequence(genome_record, cluster, 500)

# Step 2: Design sgRNAs
genome_seq = str(genome_record.seq)
sgRNAs = design_sgRNAs(cluster_seq, upstream, downstream, genome_seq, config)

# Step 3: Extract fragment
fragment = extract_fragment(genome_seq, sgRNAs["left"][0], sgRNAs["right"][0])

# Step 4: Design homology arms
arms = design_homology_arms(fragment.sequence, genome_seq, str(vector_record.seq), config)

# Step 5: Design tailed primers
primers = design_tailed_primers(arms["left"], arms["right"], str(vector_record.seq), 2300, 2306, config)

# Step 6: Simulate assembly
assembly = simulate_pydna_assembly(fragment.sequence, str(vector_record.seq),
                                    arms["left"], arms["right"], 2300, 2306)

# Step 7: Generate report
generate_report(sgRNAs, fragment.sequence, arms, primers, assembly,
                cluster, config, {"length": 5000, "gc_percent": 70.0, "valid": True},
                "results/")
```

## Output Files

The pipeline generates the following files in the output directory:

| File | Description |
|------|-------------|
| `report.md` | Markdown report with all design results |
| `report.html` | HTML report with embedded SVG cluster map |
| `cluster_map.svg` | Schematic map of the BGC with cut sites and arms |
| `primers.csv` | Primer sequences ready for synthesis ordering |

## Pipeline Steps

### Step 1: Ingestion and Annotation
- Reads BGC FASTA or GenBank file
- Validates sequence integrity (length, invalid characters, GC%)
- Identifies cluster coordinates (auto-parsed from GenBank or user-specified)

### Step 2: sgRNA Design (SpCas9 NGG)
- Scans ±500 bp (configurable) flanking regions for all NGG PAM sites
- For each candidate: calculates cut position (-3 bp from PAM), GC%, poly-T flag, and
  genomic specificity
- **BLAST+** is used if available; otherwise, a k-mer + Hamming distance fallback
  is used (flagged in report)
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
- Uses **primer3-py** for optimization:
  - Tm of annealing region (target 58-62°C)
  - Hairpin and dimer checks (primer3 thermodynamics)
  - GC clamp at the 3' end
- Cross-dimer check between forward and reverse primers

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
4. The capture vector is linearized by a **unique restriction enzyme** (default EcoRI),
   and the primers anneal to sequences flanking this cut site.
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

# Run end-to-end on synthetic data
pipeline-tar-crispr run \
    --bgc test_data/synthetic_bgc.fasta \
    --vector test_data/synthetic_vector.fasta \
    --genome test_data/synthetic_genome.fasta \
    --genbank test_data/synthetic_bgc.gbk \
    --output test_output/ \
    --verbose
```

## Fallback Behavior (When External Tools Are Unavailable)

| Feature | Primary Tool | Fallback |
|---------|-------------|----------|
| Specificity check | BLAST+ (blastn) | k-mer + Hamming distance search |
| Secondary structure | ViennaRNA (RNAfold) | Inverted-repeat heuristic |

When fallbacks are used, a warning is logged and the `blast_available`/`rnafold_available` flags are set to `False` in the report.

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
