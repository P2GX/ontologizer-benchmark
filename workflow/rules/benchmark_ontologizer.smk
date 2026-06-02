RECALL_LEVELS    = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
PRECISION_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
REPEATS          = list(range(10))
FIXED_PRECISION  = 0.5
FIXED_RECALL     = 0.5
METHODS = ["freq", "bayes"]

configfile: "config/config.yaml"

rule gen_population_genes:
    input:
        go = config["go_obo"],
        gaf = config["gaf"],
    output:
        "results/synthetic/population_genes.txt",
    script:
        "../scripts/synthetic/gen_population_genes.py"

rule gen_study_genes_recall:
    input:
        go = config["go_obo"],
        gaf = config["gaf"],
    output:
        "results/synthetic/recall/recall_{recall}.tsv",
    params:
        precision = FIXED_PRECISION,
    script:
        "../scripts/synthetic/gen_recall_study_genes.py"

rule gen_study_genes_precision:
    input:
        go = config["go_obo"],
        gaf = config["gaf"],
    output:
        "results/synthetic/precision/precision_{precision}.tsv",
    params:
        recall = FIXED_RECALL,
    script:
        "../scripts/synthetic/gen_precision_study_genes.py"

rule run_freq_benchmark:
    input:
        binary = config["ontologizer"],
        go = config["go_json"],
        gaf = config["gaf"],
        population = "results/synthetic/population_genes.txt",
        study = "results/synthetic/{stat}/{stat}_{value}.tsv",
    output:
        "results/synthetic/{stat}/results/ora_freq_{stat}_{value}_repeat_{repeat}.csv"
    script:
        "../scripts/synthetic/run_freq_benchmark.py"

rule run_bayesian_benchmark:
    input:
        binary = config["ontologizer"],
        go = config["go_json"],
        gaf = config["gaf"],
        population = "results/synthetic/population_genes.txt",
        study = "results/synthetic/{stat}/{stat}_{value}.tsv",
    output:
        "results/synthetic/{stat}/results/ora_bayes_{stat}_{value}_repeat_{repeat}.csv"
    script:
        "../scripts/synthetic/run_bayes_benchmark.py"

rule plot_benchmark:
    input:
        recall_genes    = expand("results/synthetic/recall/recall_{recall}.tsv",
            recall=RECALL_LEVELS),
        precision_genes = expand("results/synthetic/precision/precision_{precision}.tsv",
            precision=PRECISION_LEVELS),
        freq_recall     = expand("results/synthetic/recall/results/ora_freq_recall_{recall}_repeat_{repeat}.csv",
            recall=RECALL_LEVELS, repeat=REPEATS),
        bayes_recall    = expand("results/synthetic/recall/results/ora_bayes_recall_{recall}_repeat_{repeat}.csv",
            recall=RECALL_LEVELS, repeat=REPEATS),
        freq_precision  = expand("results/synthetic/precision/results/ora_freq_precision_{precision}_repeat_{repeat}.csv",
            precision=PRECISION_LEVELS, repeat=REPEATS),
        bayes_precision = expand("results/synthetic/precision/results/ora_bayes_precision_{precision}_repeat_{repeat}.csv",
            precision=PRECISION_LEVELS, repeat=REPEATS),
    output:
        png = "results/benchmark_ontologizer.png",
        pdf = "results/benchmark_ontologizer.pdf",
    params:
        recalls    = RECALL_LEVELS,
        precisions = PRECISION_LEVELS,
        n_repeats  = len(REPEATS),
    script:
        "../scripts/plot/plot_ontologizer.py"