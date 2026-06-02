# Ontologizer Benchmark

Benchmarks **Ontologizer's** GO (Gene Ontology) over-representation analysis (ORA):

- **Frequentist (TfT)** — Term-for-Term Fisher's exact test, Bonferroni-corrected. Terms called significant at *p* ≤ 0.05.
- **Bayesian (MGSA)** — Model-based Gene-set Analysis. Terms called significant at posterior probability ≥ 0.50.

The benchmark is **synthetic**: study gene sets are generated from real GO/GAF annotations with controlled signal, 
so the ground-truth of active terms is known and detection precision/recall can be measured exactly.

## Prerequisites

- Conda
- A pre-built Ontologizer binary, compiled separately from the `ontologizer`
  Rust project (https://github.com/P2GX/ontologizer). This repository only invokes it via `subprocess`; it does not
  build or ship it.

## Setup

**1. Create and activate the Conda environment:**

```bash
conda env create -f workflow/envs/environment.yaml
conda activate ontologizer-benchmark
```

**2. Point `config/config.yaml` at your Ontologizer binary:**

```yaml
go_json:     resources/go-basic.json   # auto-downloaded
go_obo:      resources/go-basic.obo    # auto-downloaded
gaf:         resources/goa_human.gaf   # auto-downloaded
ontologizer: /path/to/ontologizer      # <-- compiled Rust binary (absolute path)
```

The GO files (`go-basic.json`, `go-basic.obo`, `goa_human.gaf`) are downloaded
automatically by the `download_go` rule the first time the pipeline runs.

## Running

```bash
# Confirm the workflow graph resolves without running anything
snakemake -n

# Run the full pipeline
snakemake --cores <N> --use-conda
```

This produces two files in `results/`:

| File | Description |
|---|---|
| `benchmark_ontologizer.png` | Precision/recall plot (raster) |
| `benchmark_ontologizer.pdf` | Precision/recall plot (vector) |

## Synthetic Benchmark Design

Each study set contains roughly 500 signal genes drawn from randomly chosen GO terms, plus noise genes, and is characterized by two parameters:

- **ρ** — the fraction of each active term's annotated genes included in the study set.
- **η** — the fraction of study-set genes annotated to an active term (the remainder are noise).

For each (ρ, η) combination, 10 study sets are sampled. Every study set is scored by both Ontologizer modes, and the plotter compares each method's significant terms against the known ground truth, reporting term **precision** and **recall**.

## Project Structure

```
config/
  config.yaml                    # Paths to GO data and the Ontologizer binary
resources/                       # Downloaded GO/GAF data (git-ignored)
results/                         # Generated gene sets, ORA results, plots (git-ignored)
workflow/
  Snakefile                      # Entry point
  rules/
    download_go.smk              # Download GO hierarchy + human annotations
    benchmark_ontologizer.smk    # Synthetic benchmark rules
  scripts/
    synthetic/
      common.py                  # term_gene_map / study-set generation / extraction
      gen_population_genes.py     # Background gene universe
      gen_recall_study_genes.py   # Study sets for the recall (ρ) sweep
      gen_precision_study_genes.py# Study sets for the precision (η) sweep
      run_freq_benchmark.py       # Frequentist (FET) ORA via Ontologizer
      run_bayes_benchmark.py      # Bayesian (MGSA) ORA via Ontologizer
    plot/
      plot_ontologizer.py         # Build the 2x2 benchmark figure
    others/
      gen_population_genes_organism.py  # Non-human helpers (not in the Snakefile)
      gen_study_genes_organism.py
  envs/
    environment.yaml             # Conda environment definition
```