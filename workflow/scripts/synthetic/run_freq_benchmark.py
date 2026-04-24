import json
import os
import subprocess
import tempfile
from typing import Dict, Any

from common import extract_study_genes

if "snakemake" not in dir():
    class snakemake:
        class input:
            binary     = "/home/lukas/RustroverProjects/ontologizer/target/release/ontologizer"
            go         = "/home/lukas/RustroverProjects/ontologizer/data/go-basic.json"
            gaf        = "/home/lukas/RustroverProjects/ontologizer/data/goa_human.gaf"
            study      = "../../results/synthetic/recall/recall_0.2.tsv"
            population = "../../results/synthetic/population_genes.txt"
        class wildcards:
            fmethod = "ParentIntersection"
            stat = "recall"
            value = "0.2"
            repeat = "0"
        output = [f"../../results/synthetic/recall/results/ora_freq_ParentIntersection_recall_0.2_repeat_{i}.csv" for i in range(10)]
    os.makedirs("../../results/synthetic/recall/results", exist_ok=True)

binary = snakemake.input.binary
go = snakemake.input.go
gaf = snakemake.input.gaf
study_tsv = snakemake.input.study
population = snakemake.input.population

statistic  = snakemake.wildcards.stat
value = float(snakemake.wildcards.value)
repeat = int(snakemake.wildcards.repeat)

out_path = snakemake.output[0]

# Extract study genes for this repeat to a temporary file
with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
    study_path = tmp.name
extract_study_genes(study_tsv, repeat=repeat, out_path=study_path)

# Write ontologizer config to a temporary file
ora_config: Dict[str, Any] = {
    "study_file": study_path,
    "pop_file":   population,
    "go_file":    go,
    "goa_file":   gaf,
    "out_file":   out_path,
    "method": {
        "method":     "Frequentist",
        "background": "Standard",
        "correction": "Bonferroni",
    },
}
with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tmp:
    json.dump(ora_config, tmp, indent=4)
    config_path = tmp.name

subprocess.run([binary, config_path])

os.remove(study_path)
os.remove(config_path)