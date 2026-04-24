DATASETS = ["GSE1297", "GSE14762", "GSE15471" ,"GSE16515", "GSE18842", "GSE19188", "GSE19728", "GSE20153", "GSE20291",
            "GSE21354", "GSE3467", "GSE3585", "GSE3678", "GSE4107", "GSE5281_EC", "GSE5281_HIP", "GSE5281_VCX",
            "GSE6956AA", "GSE6956C", "GSE781", "GSE8671", "GSE8762", "GSE9348", "GSE9476", "GSE1145", "GSE11906",
            "GSE14924_CD4", "GSE14924_CD8", "GSE16759", "GSE19420", "GSE20164", "GSE22780", "GSE23878", "GSE24739_G0",
            "GSE24739_G1", "GSE30153", "GSE32676", "GSE38666_epithelia", "GSE38666_stroma", "GSE4183", "GSE42057",
            "GSE7305"]

rule download_geo2kegg:
    output:
        expand("results/geo2kegg/{set}/{set}_study.txt", set=DATASETS),
        expand("results/geo2kegg/{set}/{set}_population.txt", set=DATASETS)
    script:
        "../scripts/geo2kegg/download_geo2kegg.R"

rule build_symbol_tab:
    input:
        gaf = config["gaf"]
    output:
        tab = "resources/goa_human_symbol.tab"
    script:
        "../scripts/geo2kegg/build_symbol_tab.py"

rule run_goa:
    input:
        study = "results/geo2kegg/{set}/{set}_study.txt",
        population = "results/geo2kegg/{set}/{set}_population.txt",
        tab = "resources/goa_human_symbol.tab"
    output:
        "results/geo2kegg/{set}/{set}_freq_goatools_ora.csv"
    script:
        "../scripts/geo2kegg/run_freq_goatools.py"

rule run_freq_ontologizer:
    input:
        study = "results/geo2kegg/{set}/{set}_study.txt",
        population = "results/geo2kegg/{set}/{set}_population.txt"
    output:
        "results/geo2kegg/{set}/{set}_freq_ontologizer_ora.csv"
    script:
        "../scripts/geo2kegg/run_freq_ontologizer.py"

rule run_bayes_ontologizer:
    input:
        study = "results/geo2kegg/{set}/{set}_study.txt",
        population = "results/geo2kegg/{set}/{set}_population.txt"
    output:
        "results/geo2kegg/{set}/{set}_bayes_ontologizer_ora.csv"
    script:
        "../scripts/geo2kegg/run_bayes_ontologizer.py"

rule plot_geo2kegg_benchmark:
    input:
        goa = expand("results/geo2kegg/{set}/{set}_freq_goatools_ora.csv", set=DATASETS),
        freq = expand("results/geo2kegg/{set}/{set}_freq_ontologizer_ora.csv", set=DATASETS),
        bayes = expand("results/geo2kegg/{set}/{set}_bayes_ontologizer_ora.csv", set=DATASETS),
    output:
        png = "results/benchmark_geo2kegg.png",
        pdf = "results/benchmark_geo2kegg.pdf",
    script:
        "../scripts/plot/plot_geo2kegg.py"