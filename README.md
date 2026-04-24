# Ontologizer Benchmark

Benchmarks and compares GO (Gene Ontology) enrichment analysis methods:

- **Ontologizer** — a Rust-based ORA tool supporting frequentist (Bonferroni-corrected) and Bayesian scoring
- **GOATools** — a Python-based GO enrichment library

Two evaluation tracks are included:

1. **Synthetic benchmark** — controlled recall and noise experiments with known ground truth
2. **geo2kegg benchmark** — 42 real-world GEO datasets with KEGG pathway annotations

## Prerequisites

- Conda
- A built Ontologizer binary (from the `ontologizer` Rust project)
- GO files: `go-basic.obo`, `go-basic.json`, `goa_human.gaf`

## Setup

**1. Create the Conda environment:**

```bash
conda env create -f workflow/envs/environment.yaml
conda activate <env-name>
```

**2. Configure paths in `config/config.yaml`:**

```yaml
go_json:      /path/to/go-basic.json
go_obo:       /path/to/go-basic.obo
gaf:          /path/to/goa_human.gaf
ontologizer:  /path/to/ontologizer   # compiled Rust binary
```

**3. Download the geo2kegg datasets** (requires R + BioConductor):

```bash
snakemake --cores 1 --use-conda download_geo2kegg
```

## Running

Run the full pipeline:

```bash
snakemake --cores <N> --use-conda
```

This produces four output files in `results/`:

| File | Description |
|---|---|
| `benchmark_ontologizer.png/pdf` | Precision/recall curves for the synthetic benchmark |
| `benchmark_geo2kegg.png/pdf` | Method comparison on real GEO datasets |

## Project Structure

```
config/
  config.yaml                  # Paths to GO data and Ontologizer binary
resources/
  go-basic.obo                 # GO hierarchy (OBO format)
  go-basic.json                # GO hierarchy (JSON format)
  goa_human.gaf                # Human gene-GO annotations
  goa_human_symbol.tab         # Gene symbol lookup table (generated)
results/
  synthetic/                   # Generated synthetic gene sets and ORA results
  geo2kegg/                    # GEO datasets and ORA results
workflow/
  Snakefile
  rules/
    benchmark_ontologizer.smk  # Synthetic benchmark rules
    benchmark_geo2kegg.smk     # geo2kegg benchmark rules
  scripts/
    common.py                  # Shared utilities
    gen_population_genes.py    # Background gene universe
    gen_recall_study_genes.py  # Synthetic sets varying signal recall
    gen_noise_study_genes.py   # Synthetic sets varying noise level
    run_freq_benchmark.py      # Frequentist ORA on synthetic data
    run_bayes_benchmark.py     # Bayesian ORA on synthetic data
    run_freq_ontologizer.py    # Frequentist ORA on GEO data
    run_bayes_ontologizer.py   # Bayesian ORA on GEO data
    run_freq_goatools.py       # GOATools ORA on GEO data
    plot_ontologizer.py        # Plot synthetic benchmark results
    plot_geo2kegg.py           # Plot geo2kegg benchmark results
    download_geo2kegg.R        # Download GEO datasets (R)
    build_symbol_tab.py        # Build gene symbol lookup table
  envs/
    environment.yaml           # Conda environment definition
```

## Synthetic Benchmark Design

Synthetic gene sets are generated with controlled parameters:

- **Recall experiment** — noise fixed at 0.4, signal recall varied from 0.1 to 1.0 in steps of 0.1
- **Noise experiment** — recall fixed at 0.4, noise varied from 0.0 to 2.0× signal size in steps of 0.2
- 10 replicates per configuration, 10 recall levels × 11 noise levels = ~400 total ORA runs

Each run is evaluated by comparing predicted significant terms against the known ground-truth terms using detection recall and precision.
