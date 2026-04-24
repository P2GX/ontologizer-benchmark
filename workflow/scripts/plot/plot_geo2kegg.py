import os.path
from typing import Set, Tuple, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


if "snakemake" not in dir():
    _ROOT = "/home/lukas/PycharmProjects/Ontologizer-Benchmar"
    _DATASETS = ["GSE1297", "GSE14762"]
    _GEO = f"{_ROOT}/resources/geo2kegg"

    class snakemake:
        class input:
            goa   = [f"{_GEO}/{s}/{s}_freq_goatools_ora.csv"    for s in _DATASETS]
            freq  = [f"{_GEO}/{s}/{s}_freq_ontologizer_ora.csv" for s in _DATASETS]
            bayes = [f"{_GEO}/{s}/{s}_bayes_ontologizer_ora.csv" for s in _DATASETS]
        class output:
            png = f"{_ROOT}/results/benchmark_geo2kegg.png"
            pdf = f"{_ROOT}/results/benchmark_geo2kegg.pdf"


def calculate_hierarchical_bars(set_a: Set, set_b: Set, set_c: Set):
    """
    Computes the cardinalities of five mutually disjoint topological regions.

    Mapping Context:
    A = Bayesian
    B = Frequentist
    C = Goa
    """

    # Pre-compute structural boundaries
    union_bc = set_b | set_c
    union_ac = set_a | set_c
    union_ab = set_a | set_b

    # 1. Elements shared by all three sets: A ∩ B ∩ C
    intersection_all = len(set_a & set_b & set_c)

    # 2. Elements exclusively in A: A \ (B ∪ C)
    exclusive_A = len(set_a - union_bc)

    # 3. Elements shared by B and C, but strictly missing from A: (B ∩ C) \ A
    shared_B_C_only = len((set_b & set_c) - set_a)

    # 4. Elements exclusively in B: B \ (A ∪ C)
    exclusive_B = len(set_b - union_ac)

    # 5. Elements exclusively in C: C \ (A ∪ B)
    exclusive_C = len(set_c - union_ab)

    return [intersection_all, exclusive_A, shared_B_C_only, exclusive_B, exclusive_C]


def plot_hierarchical_bars(
    hierarchy: List[Tuple[int, int, int, int, int]],
    labels: List[str],
    out_png: str,
    out_pdf: str,
) -> None:
    """
    Generates a stacked bar chart for an arbitrary sequence of 5-region sets.

    Args:
        hierarchy: A list of 5-element tuples corresponding to:
                   [intersection_all, exclusive_A, shared_B_C_only, exclusive_B, exclusive_C]
        labels:    Dataset name for each bar.
        out_png:   File path for the PNG output.
        out_pdf:   File path for the PDF output.
    """
    # Convert sequence to a 2D numpy array for vectorized column operations
    data = np.array(hierarchy)

    if data.shape[1] != 5:
        raise ValueError(f"Expected 5 regions per item, received {data.shape[1]}")

    num_bars = data.shape[0]
    x_positions = np.arange(num_bars)

    # Extract structural regions into independent 1D arrays
    intersection_all = data[:, 0]
    exclusive_A = data[:, 1]
    shared_B_C_only = data[:, 2]
    exclusive_B = data[:, 3]
    exclusive_C = data[:, 4]

    # Calculate cumulative offsets for the stacked bars
    bottom_exclusive_A = intersection_all
    bottom_shared_BC = bottom_exclusive_A + exclusive_A
    bottom_exclusive_B = bottom_shared_BC + shared_B_C_only
    bottom_exclusive_C = bottom_exclusive_B + exclusive_B

    fig, ax = plt.subplots(figsize=(8, 6))
    width = 0.6

    # Define a distinct 5-color categorical palette
    colors = ["#E69F00", "#298C8C", "#0072B2", "#298C8C", "#a00000"]

    # Plot each layer using the pre-calculated bottoms
    ax.bar(x_positions, intersection_all, width,
           label="Shared All", color=colors[0], edgecolor="white", linewidth=1)
    ax.bar(x_positions, exclusive_A, width, bottom=bottom_exclusive_A,
           label="Bayes not in Frequentist", color=colors[1], edgecolor="white", linewidth=1)
    ax.bar(x_positions, shared_B_C_only, width, bottom=bottom_shared_BC,
           label="Frequentist not in Bayes", color=colors[2], edgecolor="white", linewidth=1)
    ax.bar(x_positions, exclusive_B, width, bottom=bottom_exclusive_B,
           label="Frequentist not in Bayes and Ontologizer only", color=colors[3], edgecolor="white", linewidth=1)
    ax.bar(x_positions, exclusive_C, width, bottom=bottom_exclusive_C,
           label="Frequentist not in Bayes and GOA only", color=colors[4], edgecolor="white", linewidth=1)

    # Formatting
    ax.set_ylabel("Significant terms")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()


def goa_significance(goa_results) -> Set:
    """Filter GOA results to only include significant enriched terms."""
    x = goa_results[(goa_results["enrichment"] == "e") & (goa_results["p_bonferroni"] <= 0.05)]
    return set(x["term_id"])


def onto_freq_significance(onto_freq_results) -> Set:
    """Filter Ontologizer frequentist results to only include significant terms."""
    x = onto_freq_results[onto_freq_results["Score"] <= 0.05]
    return set(x["Id"])


def onto_bayesian_significance(onto_bayes_results) -> Set:
    """Filter Ontologizer Bayesian results to only include significant terms."""
    x = onto_bayes_results[onto_bayes_results["Score"] >= 0.50]
    return set(x["Id"])


overlap = dict()
for goa_path, freq_path, bayes_path in zip(
    snakemake.input.goa, snakemake.input.freq, snakemake.input.bayes
):
    name = os.path.basename(os.path.dirname(goa_path))

    goa_results        = pd.read_csv(goa_path)
    onto_freq_results  = pd.read_csv(freq_path,  skiprows=3)
    onto_bayes_results = pd.read_csv(bayes_path, skiprows=4 )

    goa_sig        = goa_significance(goa_results)
    onto_freq_sig  = onto_freq_significance(onto_freq_results)
    onto_bayes_sig = onto_bayesian_significance(onto_bayes_results)
    print(name, len(goa_sig), len(onto_freq_sig), len(onto_bayes_sig))
    overlap[name] = calculate_hierarchical_bars(onto_bayes_sig, onto_freq_sig, goa_sig)

plot_hierarchical_bars(
    list(overlap.values()),
    list(overlap.keys()),
    snakemake.output.png,
    snakemake.output.pdf,
)