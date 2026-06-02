import json
import os
import subprocess
import tempfile
from typing import Dict, Any

from common import extract_study_genes

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