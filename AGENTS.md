# AGENTS.md

Project-specific guidance for working with this codebase.

## Testing Commands

All tests:
```bash
cd C:\kilo.tar-crispr-pipeline
python -m pytest tests/ -v
```

Run specific module tests:
```bash
python -m pytest tests/test_pam_finder.py tests/test_homology_arms.py -v
```

End-to-end CLI test (also included as a pytest test):
```bash
python -m tar_crispr.cli run --bgc test_data/synthetic_bgc.fasta --vector test_data/synthetic_vector.fasta --genome test_data/synthetic_genome.fasta --genbank test_data/synthetic_bgc.gbk --output test_output/ --verbose --no-blast
```

## Architecture

- **Language:** Python 3.10+
- **Package:** `tar_crispr` under `src/`
- **CLI:** Typer-based, entry point `pipeline-tar-crispr`
- **Key libraries:** Biopython, primer3-py, pydna, matplotlib, typer

## Module Overview

| Module | Step | Responsibility |
|--------|------|----------------|
| `config.py` | — | Dataclasses: PipelineConfig, ClusterInfo, PAMCandidate, HomologyArm, TailedPrimer, AssemblyResult |
| `sequence_io.py` | 1 | FASTA/GenBank parsing, sequence validation, cluster coordinate extraction |
| `pam_finder.py` | 2 | PAM scanning, sgRNA scoring/ranking, specificity (BLAST or k-mer fallback) |
| `fragment_ends.py` | 3 | Cas9 blunt-cut position calculation, fragment extraction |
| `homology_arms.py` | 4 | Arm extraction, validation (GC, uniqueness, secondary structure, RE sites), shifting |
| `primer_design.py` | 5 | Tailed primer construction, primer3 optimization (Tm, hairpin, dimer, GC clamp) |
| `assembly_sim.py` | 6 | Assembly simulation (pydna + string fallback), junction/spurious overlap checks |
| `report_generator.py` | 7 | Markdown report, SVG cluster map, HTML report, CSV primer export |
| `cli.py` | — | Typer CLI orchestrating all steps |

## Key Conventions

1. **Coordinates are 0-based** throughout the codebase (Python/BioPython convention).
2. **All sequences are uppercased** before processing (handled in each module).
3. **BLAST/RNAfold are optional** — the pipeline auto-detects availability and falls back gracefully.
4. **No hardcoded paths** — everything is CLI-configurable or uses sensible defaults.
5. **Cut positions** are always relative to the full genomic sequence (offsets applied in `design_sgRNAs`).

## Test Files

| File | Coverage |
|------|----------|
| `tests/conftest.py` | Shared fixtures (synthetic sequences, temp file creators) |
| `tests/test_sequence_io.py` | FASTA/GenBank parsing, validation, cluster extraction |
| `tests/test_pam_finder.py` | PAM detection, GC calc, specificity, ranking, sgRNA design |
| `tests/test_fragment_ends.py` | Cut position calculation, fragment extraction |
| `tests/test_homology_arms.py` | Arm extraction, GC validation, uniqueness, secondary structure, shifting |
| `tests/test_primer_design.py` | Tm, GC, hairpin, dimer, tailed primer construction, CSV export |
| `tests/test_assembly_sim.py` | Junction checks, spurious overlaps, assembly simulation |
| `tests/test_report_generator.py` | Table formatting, SVG generation, full report generation |
| `tests/test_e2e.py` | End-to-end CLI runs on synthetic test data (4 tests) |
