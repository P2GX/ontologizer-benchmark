RECALL_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
NOISE_LEVELS  = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
REPEATS       = list(range(10))
FIXED_NOISE   = 0.4
FIXED_RECALL  = 0.4
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
        noise = FIXED_NOISE,
    script:
        "../scripts/synthetic/gen_recall_study_genes.py"

rule gen_study_genes_noise:
    input:
        go = config["go_obo"],
        gaf = config["gaf"],
    output:
        "results/synthetic/noise/noise_{noise}.tsv",
    params:
        recall = FIXED_RECALL,
    script:
        "../scripts/synthetic/gen_noise_study_genes.py"

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

rule plot_ontologizer_benchmark:
    input:
        recall_genes = expand("results/synthetic/recall/recall_{recall}.tsv",
            recall=RECALL_LEVELS),
        noise_genes  = expand("results/synthetic/noise/noise_{noise}.tsv",
            noise=NOISE_LEVELS),
        freq_recall  = expand("results/synthetic/recall/results/ora_freq_recall_{recall}_repeat_{repeat}.csv",
            recall=RECALL_LEVELS, repeat=REPEATS),
        bayes_recall = expand("results/synthetic/recall/results/ora_bayes_recall_{recall}_repeat_{repeat}.csv",
            recall=RECALL_LEVELS, repeat=REPEATS),
        freq_noise   = expand("results/synthetic/noise/results/ora_freq_noise_{noise}_repeat_{repeat}.csv",
            noise=NOISE_LEVELS, repeat=REPEATS),
        bayes_noise  = expand("results/synthetic/noise/results/ora_bayes_noise_{noise}_repeat_{repeat}.csv",
            noise=NOISE_LEVELS, repeat=REPEATS),
    output:
        png = "results/benchmark_ontologizer.png",
        pdf = "results/benchmark_ontologizer.pdf",
    params:
        recalls   = RECALL_LEVELS,
        noises    = NOISE_LEVELS,
        n_repeats = len(REPEATS),
    script:
        "../scripts/plot/plot_ontologizer.py"

