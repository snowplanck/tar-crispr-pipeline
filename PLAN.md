# Implementation Plan: TAR-CRISPR Pipeline

## Environment Status
- Python 3.14.7
- Available: biopython, primer3-py, pydna, typer, pytest, matplotlib
- **Not available:** BLAST+ (`blastn`), ViennaRNA (`RNAfold`)
- Therefore: fallbacks for specificity (k-mer/regex) and secondary structure (inverted-repeat heuristic) are mandatory

## Project Structure
```
tar_crispr_pipeline/
├── pyproject.toml
├── README.md
├── AGENTS.md
├── src/
│   └── tar_crispr/
│       ├── __init__.py
│       ├── cli.py              # Typer CLI entry point
│       ├── config.py           # Dataclasses / config defaults
│       ├── sequence_io.py      # Step 1: Ingestion, validation, cluster coords
│       ├── pam_finder.py       # Step 2: PAM/sgRNA discovery, scoring, ranking
│       ├── fragment_ends.py    # Step 3: Cut calculation, fragment extraction
│       ├── homology_arms.py    # Step 4: Arm extraction, validation, shifting
│       ├── primer_design.py    # Step 5: Tailed primers via primer3-py
│       ├── assembly_sim.py     # Step 6: In silico assembly simulation
│       └── report_generator.py # Step 7: Markdown + SVG report, CSV export
├── tests/
│   ├── conftest.py             # Pytest fixtures (synthetic data)
│   ├── test_pam_finder.py      # PAM detection, cut position, ranking
│   ├── test_homology_arms.py   # Arm extraction, shift logic, validation
│   ├── test_primer_design.py   # Tm validation, primer generation
│   ├── test_assembly_sim.py    # Assembly simulation correctness
│   └── test_sequence_io.py     # FASTA/GenBank parsing, validation
└── test_data/
    ├── synthetic_bgc.fasta     # Fictitious BGC (~3kb)
    ├── synthetic_vector.fsa      # Capture vector (~5kb)
    ├── synthetic_genome.fasta    # Genome containing BGC for specificity
    └── synthetic_bgc.gbk         # GenBank with gene annotations
```

## Module-Level Implementation Details

### 1. config.py
- `PipelineConfig` dataclass with all defaults: PAM window 500bp, arm length 50bp,
  Tm range 58-62C, max shift 100bp, top N=5, mismatch tolerance 3, etc.
- All parameters overridable via CLI flags.

### 2. sequence_io.py (Step 1)
- `read_sequence(path)` — auto-detects FASTA or GenBank via Biopython.
- `validate_sequence(seq)` — checks for invalid chars, reports length, GC%.
- `extract_cluster_bounds(record, start=None, end=None)` — if start/end not given,
  parse from GenBank features (e.g. `feature.type == "cluster"` or CDS span).
- Returns `SeqRecord` and `ClusterInfo` dataclass with coordinates.

### 3. pam_finder.py (Step 2)
- `find_pam_sites(seq, pam="NGG", strand="+")` — scans both strands, returns list
  of `PAMCandidate` dataclass: position, strand, protospacer_seq, pam_seq, cut_pos.
- `calculate_cut_position(pam_pos, strand)` — SpCas9 cuts at -3 from NGG (blunt).
- `calculate_gc(seq)` — simple GC%.
- `has_polyt(seq, k=4)` — flags 4+ consecutive T's.
- `check_specificity(protospacer, genome_seq, max_mismatch=3)` — BLAST if available,
  otherwise k-mer + Hamming distance fallback.
- `rank_sgRNAs(candidates, genome_seq)` — scores and returns top N per end.
- Public API: `design_sgRNAs(cluster_seq, left_flank, right_flank, genome_seq, config)`.

### 4. fragment_ends.py (Step 3)
- `compute_fragment(sg_left, sg_right, bgc_seq)` — determines the cut fragment
  from left cut position to right cut position across the BGC.
- Handles strand orientation: sgRNA on + strand vs - strand changes cut coordinates.

### 5. homology_arms.py (Step 4)
- `extract_homology_arm(fragment, end, length=50)` — extracts arm from cut end.
- `validate_arm(arm, genome_seq, vector_seq, avoid_enzymes=None)` — checks
  uniqueness (k-mer fallback), GC range, secondary structure (inverted-repeat),
  restriction sites.
- `find_valid_arm(fragment, end, config, genome_seq, vector_seq)` — tries default
  length, then shifts in 5bp increments up to max_shift.
- Public API: `design_homology_arms(fragment, genome_seq, vector_seq, config)`.

### 6. primer_design.py (Step 5)
- `design_tailed_primer(arm, vector_seq, vector_cut_pos, end, config)` — builds
  primer = arm + annealing_region; uses primer3-py to optimize annealing region.
- `validate_primer(primer_seq)` — Tm, hairpin, dimer checks (primer3 built-in +
  manual GC clamp check).
- Returns `TailedPrimer` dataclass with full seq, tail, annealing region, metrics.
- Exports to CSV via `export_primers(primers, path)`.

### 7. assembly_sim.py (Step 6)
- `simulate_assembly(fragment, vector_linearized, left_arm, right_arm)` — uses pydna
  logic: simulate Gibson/HR assembly, verify junctions match, circle closes, no
  unexpected overlaps, size check.
- Uses `pydna.DNAseq` for manipulation or manual string operations as fallback.
- Returns `AssemblyResult` with success flag, final size, any issues.

### 8. report_generator.py (Step 7)
- `generate_report(design, output_dir)` — creates Markdown report.
- SVG diagram of BGC with cut sites, arms, primers marked.
- Tables for sgRNAs, homology arms, primers.
- CSV export of primers.
- Uses matplotlib for simple coordinate map (renders to SVG).

### 9. cli.py
- Typer app: `pipeline-tar-crispr --bgc ... --vector ... --output ...`
- Optional flags: `--genome`, `--genbank`, `--start`, `--end`, `--vector-cut-site`,
  `--top-n`, `--arm-length`, `--pam-window`, `--no-blast` (force fallback mode).
- Orchestrates all steps sequentially, saving intermediate results.

## Test Strategy
- `conftest.py`: generates in-memory synthetic sequences (BGC ~3kb, vector ~5kb,
  genome ~6kb containing BGC) for reuse across tests.
- `test_pam_finder.py`: known PAM sequences, cut position verification, poly-T flag,
  GC%, ranking.
- `test_homology_arms.py`: arm extraction at exact coordinates, shift logic,
  validation of uniqueness/GC/secondary structure.
- `test_primer_design.py`: primer3 Tm validation, tailed primer construction.
- `test_assembly_sim.py`: assembly of known fragment+vector, circle closure check.
- End-to-end: run full CLI on synthetic test_data files, verify report exists and
  contains all sections.

## Execution Order (test each module before integrating)
1. config.py → 2. sequence_io.py (test: parsing/validation) →
3. pam_finder.py (test: PAM, cut, scoring) →
4. fragment_ends.py (test: cut fragment extraction) →
5. homology_arms.py (test: arm extraction/validation/shift) →
6. primer_design.py (test: tailed primer/Tm) →
7. assembly_sim.py (test: assembly simulation) →
8. report_generator.py →
9. cli.py →
10. End-to-end run with synthetic data
