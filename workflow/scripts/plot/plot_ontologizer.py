import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from matplotlib.axes import Axes

# ── Publication-ready Matplotlib settings ─────────────────────────────────
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 0.75,
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})

# ── Okabe-Ito colorblind-safe palette ─────────────────────────────────────────
FREQ_COLOR = "#0072B2"  # blue
BAYES_COLOR = "#E69F00"  # orange


def load_true_terms(file, repeat):
    data = pd.read_csv(file, sep="\t")
    data_n = data[(data["repeat"] == repeat) & (data["term_id"] != "Noise")]
    return set(data_n["term_id"].unique())


def load_predicted_freq_terms(file):
    data = pd.read_csv(file, sep=",", comment="!")
    return set(data[data["Score"] <= 0.05]["Id"].unique())


def load_predicted_bayesian_terms(file):
    data = pd.read_csv(file, sep=",", comment="!")
    return set(data[data["Score"] >= 0.50]["Id"].unique())


def compute_metrics(true_terms, sig_terms):
    tp = len(true_terms & sig_terms)
    detection_recall = tp / len(true_terms) if true_terms else 0.0
    precision = tp / len(sig_terms) if sig_terms else 0.0
    return detection_recall, precision


def get_method_data(gene_files, method_files, n_repeats, load_method_func):
    """Compute per-x summary statistics from pre-resolved file lists.

    gene_files   : one TSV per x value (in order)
    method_files : n_repeats CSVs per x value, grouped by x (expand order)
    """
    mean_recalls, std_recalls = [], []
    mean_precisions, std_precision = [], []

    for i, gene_file in enumerate(gene_files):
        recalls_for_x = []
        precisions_for_x = []

        for rep in range(n_repeats):
            true_terms = load_true_terms(gene_file, rep)
            method_terms = load_method_func(method_files[i * n_repeats + rep])

            rec, prec = compute_metrics(true_terms, method_terms)
            recalls_for_x.append(rec)
            precisions_for_x.append(prec)

        mean_recalls.append(np.mean(recalls_for_x))
        std_recalls.append(np.std(recalls_for_x))
        mean_precisions.append(np.mean(precisions_for_x))
        std_precision.append(np.std(precisions_for_x))


    return {
        "recall": {"mean": mean_recalls, "std": std_recalls},
        "precision": {"mean": mean_precisions, "std": std_precision},
    }


# ── Load Data ─────────────────────────────────────────────────────────────────
n_repeats  = snakemake.params.n_repeats
recalls    = snakemake.params.recalls
precisions = snakemake.params.precisions

freq_recall_files = snakemake.input.freq_recall
bayes_recall_files = snakemake.input.bayes_recall

freq_precision_files = snakemake.input.freq_precision
bayes_precision_files = snakemake.input.bayes_precision

freq_recall = get_method_data(
    snakemake.input.recall_genes, freq_recall_files, n_repeats, load_predicted_freq_terms
)

bayes_recall = get_method_data(
    snakemake.input.recall_genes, bayes_recall_files, n_repeats, load_predicted_bayesian_terms
)

freq_precision = get_method_data(
    snakemake.input.precision_genes, freq_precision_files, n_repeats, load_predicted_freq_terms
)

bayes_precision = get_method_data(
    snakemake.input.precision_genes, bayes_precision_files, n_repeats, load_predicted_bayesian_terms
)

# ── Setup Figure ──────────────────────────────────────────────────────────
# Standard 2-column width is ~170mm (6.7 inches)
fig, axes = plt.subplots(2, 2, sharex="col", sharey=True, gridspec_kw={"wspace": 0.1, "hspace": 0.1})

ax_rr : Axes = axes[0, 0]
ax_nr : Axes = axes[0, 1]  # Top row: Recall metric
ax_rp : Axes = axes[1, 0]
ax_np : Axes = axes[1, 1]  # Bottom row: Precision metric

# ── Plot Recall Experiment (Column 1) ─────────────────────────────────────
# Recall
ax_rr.errorbar(x=recalls, y=freq_recall["recall"]["mean"], yerr=freq_recall["recall"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=FREQ_COLOR)
ax_rr.errorbar(x=recalls, y=bayes_recall["recall"]["mean"], yerr=bayes_recall["recall"]["std"], capsize=3,marker='o', linestyle="none", markersize=3, color=BAYES_COLOR)

# Precision
ax_rp.errorbar(x=recalls, y=freq_recall["precision"]["mean"], yerr=freq_recall["precision"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=FREQ_COLOR)
ax_rp.errorbar(x=recalls, y=bayes_recall["precision"]["mean"], yerr=bayes_recall["precision"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=BAYES_COLOR)

# ── Plot Precision Experiment (Column 2) ──────────────────────────────────
# Recall
ax_nr.errorbar(x=precisions, y=freq_precision["recall"]["mean"], yerr=freq_precision["recall"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=FREQ_COLOR)
ax_nr.errorbar(x=precisions, y=bayes_precision["recall"]["mean"], yerr=bayes_precision["recall"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=BAYES_COLOR)

# Precision Metric
ax_np.errorbar(x=precisions, y=freq_precision["precision"]["mean"], yerr=freq_precision["precision"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=FREQ_COLOR)
ax_np.errorbar(x=precisions, y=bayes_precision["precision"]["mean"], yerr=bayes_precision["precision"]["std"], capsize=3, marker='o', linestyle="none", markersize=3, color=BAYES_COLOR)



# ── Formatting ────────────────────────────────────────────────────────────
for ax in axes.flat:
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.6)  # Readability enhancement

# Bottom Row X-Labels
ax_rp.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax_rp.set_xlabel(r"$\rho$")

ax_np.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax_np.set_xlabel(r"$\eta$")

# Left Column Y-Labels
ax_rr.set_ylabel(r"Recall")
ax_rp.set_ylabel(r"Precision")

# ── Legend ────────────────────────────────────────────────────────────────
legend_handles = [
    mlines.Line2D(
        [], [], color=FREQ_COLOR, linestyle="-",
        label="FET (p-value ≤ 0.05)"
    ),
    mlines.Line2D(
        [], [], color=BAYES_COLOR, linestyle="-",
        label="MGSA (post. probability ≥ 0.50)"
    ),
]
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.975),
    ncols=2,
    frameon=False,
)

# ── Save Outputs ──────────────────────────────────────────────────────────
# High-impact journals prefer vector graphics (.pdf/.eps) for final typesetting
fig.savefig(snakemake.output.png, dpi=300, bbox_inches="tight")
fig.savefig(snakemake.output.pdf, bbox_inches="tight")
plt.close()