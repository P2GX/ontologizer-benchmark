import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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
    data = pd.read_csv(file, sep=",", skiprows=3)
    return set(data[data["Score"] <= 0.05]["Id"].unique())


def load_predicted_bayesian_terms(file):
    data = pd.read_csv(file, sep=",", skiprows=4)
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
    mean_recalls, low_recalls, high_recalls = [], [], []
    mean_precisions, low_precisions, high_precisions = [], [], []

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
        low_recalls.append(np.percentile(recalls_for_x, 25))
        high_recalls.append(np.percentile(recalls_for_x, 75))

        mean_precisions.append(np.mean(precisions_for_x))
        low_precisions.append(np.percentile(precisions_for_x, 25))
        high_precisions.append(np.percentile(precisions_for_x, 75))

    return {
        "recall": {"mean": mean_recalls, "low": low_recalls, "high": high_recalls},
        "precision": {
            "mean": mean_precisions,
            "low": low_precisions,
            "high": high_precisions,
        },
    }


# ── Load Data ─────────────────────────────────────────────────────────────────
n_repeats = snakemake.params.n_repeats
recalls   = snakemake.params.recalls
noises    = snakemake.params.noises

freq_recall_files = snakemake.input.freq_recall
bayes_recall_files = snakemake.input.bayes_recall

freq_noise_files = snakemake.input.freq_noise
bayes_noise_files = snakemake.input.bayes_noise

freq_recall = get_method_data(
    snakemake.input.recall_genes, freq_recall_files, n_repeats, load_predicted_freq_terms
)

bayes_recall = get_method_data(
    snakemake.input.recall_genes, bayes_recall_files, n_repeats, load_predicted_bayesian_terms
)

freq_noise = get_method_data(
    snakemake.input.noise_genes, freq_noise_files, n_repeats, load_predicted_freq_terms
)

bayes_noise = get_method_data(
    snakemake.input.noise_genes, bayes_noise_files, n_repeats, load_predicted_bayesian_terms
)

# ── Setup Figure ──────────────────────────────────────────────────────────
# Standard 2-column width is ~170mm (6.7 inches)
fig, axes = plt.subplots(2, 2, sharex="col", sharey=True, figsize=(6.7, 4.5))

ax_rr, ax_nr = axes[0, 0], axes[0, 1]  # Top row: Recall metric
ax_rp, ax_np = axes[1, 0], axes[1, 1]  # Bottom row: Precision metric

# ── Plot Recall Experiment (Column 1) ─────────────────────────────────────
# Recall Metric

ax_rr.plot(
    recalls, freq_recall["recall"]["mean"],
    color=FREQ_COLOR, marker="o", linestyle="-"
)
ax_rr.fill_between(
    recalls, freq_recall["recall"]["low"], freq_recall["recall"]["high"],
    color=FREQ_COLOR, alpha=0.15, linewidth=0
)

ax_rr.plot(
    recalls, bayes_recall["recall"]["mean"],
    color=BAYES_COLOR, marker="o", linestyle="-"
)
ax_rr.fill_between(
    recalls, bayes_recall["recall"]["low"], bayes_recall["recall"]["high"],
    color=BAYES_COLOR, alpha=0.15, linewidth=0
)

# Precision Metric
ax_rp.plot(
    recalls, freq_recall["precision"]["mean"],
    color=FREQ_COLOR, marker="o", linestyle="-"
)
ax_rp.fill_between(
    recalls, freq_recall["precision"]["low"], freq_recall["precision"]["high"],
    color=FREQ_COLOR, alpha=0.15, linewidth=0
)

ax_rp.plot(
    recalls, bayes_recall["precision"]["mean"],
    color=BAYES_COLOR, marker="o", linestyle="-"
)
ax_rp.fill_between(
    recalls, bayes_recall["precision"]["low"], bayes_recall["precision"]["high"],
    color=BAYES_COLOR, alpha=0.15, linewidth=0
)

# ── Plot Noise Experiment (Column 2) ──────────────────────────────────────
# Recall Metric
ax_nr.plot(
    noises, freq_noise["recall"]["mean"],
    color=FREQ_COLOR, marker="o", linestyle="-"
)
ax_nr.fill_between(
    noises, freq_noise["recall"]["low"], freq_noise["recall"]["high"],
    color=FREQ_COLOR, alpha=0.15, linewidth=0
)

ax_nr.plot(
    noises, bayes_noise["recall"]["mean"],
    color=BAYES_COLOR, marker="o", linestyle="-"
)
ax_nr.fill_between(
    noises, bayes_noise["recall"]["low"], bayes_noise["recall"]["high"],
    color=BAYES_COLOR, alpha=0.15, linewidth=0
)

# Precision Metric
ax_np.plot(
    noises, freq_noise["precision"]["mean"],
    color=FREQ_COLOR, marker="o", linestyle="-"
)
ax_np.fill_between(
    noises, freq_noise["precision"]["low"], freq_noise["precision"]["high"],
    color=FREQ_COLOR, alpha=0.15, linewidth=0
)

ax_np.plot(
    noises, bayes_noise["precision"]["mean"],
    color=BAYES_COLOR, marker="o", linestyle="-"
)
ax_np.fill_between(
    noises, bayes_noise["precision"]["low"], bayes_noise["precision"]["high"],
    color=BAYES_COLOR, alpha=0.15, linewidth=0
)

# ── Formatting ────────────────────────────────────────────────────────────
for ax in axes.flat:
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.6)  # Readability enhancement

# Bottom Row X-Labels
ax_rp.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax_rp.set_xlabel("Gene recall")

ax_np.set_xticks([0.0, 0.4, 0.8, 1.2, 1.6, 2.0])
ax_np.set_xlabel("Noise level")

# Left Column Y-Labels
ax_rr.set_ylabel("Recall")
ax_rp.set_ylabel("Precision")

# Panel Labels (A, B)
ax_rr.text(
    -0.18, 1.05, "A", transform=ax_rr.transAxes,
    fontweight="bold", va="bottom", ha="left"
)
ax_nr.text(
    -0.05, 1.05, "B", transform=ax_nr.transAxes,
    fontweight="bold", va="bottom", ha="left"
)

fig.subplots_adjust(hspace=0.1, wspace=0.15)

# ── Legend ────────────────────────────────────────────────────────────────
legend_handles = [
    mlines.Line2D(
        [], [], color=FREQ_COLOR, linestyle="-",
        label="Frequentist (p ≤ 0.05)"
    ),
    mlines.Line2D(
        [], [], color=BAYES_COLOR, linestyle="-",
        label="Bayesian (posterior ≥ 0.50)"
    ),
]
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.04),
    ncols=2,
    frameon=False,
)

# ── Save Outputs ──────────────────────────────────────────────────────────
# High-impact journals prefer vector graphics (.pdf/.eps) for final typesetting
fig.savefig(snakemake.output.png, dpi=300, bbox_inches="tight")
fig.savefig(snakemake.output.pdf, bbox_inches="tight")
plt.close()