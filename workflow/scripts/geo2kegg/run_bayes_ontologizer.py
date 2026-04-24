import json
import os
import subprocess
import tempfile
from typing import Any, Dict

if "snakemake" not in dir():
    _ROOT = "/home/lukas/PycharmProjects/Ontologizer-Utils"
    _SET  = "GSE1297"
    _GEO  = f"{_ROOT}/resources/geo2kegg/{_SET}"

    class snakemake:
        config = {
            "ontologizer": "/home/lukas/RustroverProjects/ontologizer/target/release/ontologizer",
            "go_json":     f"{_ROOT}/resources/go-basic.json",
            "gaf":         f"{_ROOT}/resources/goa_human.gaf",
        }
        class input:
            study      = f"{_GEO}/{_SET}_study.txt"
            population = f"{_GEO}/{_SET}_population.txt"
        output = [f"{_GEO}/{_SET}_bayes_ontologizer_ora.csv"]

binary  = snakemake.config["ontologizer"]
go_json = snakemake.config["go_json"]
gaf     = snakemake.config["gaf"]

out_path     = snakemake.output[0]
dataset_dir  = os.path.dirname(out_path)
dataset_name = os.path.basename(dataset_dir)

ora_config: Dict[str, Any] = {
    "study_file": snakemake.input.study,
    "pop_file":   snakemake.input.population,
    "go_file":    go_json,
    "goa_file":   gaf,
    "out_file":   out_path,
    "method": {
        "method": "Bayesian",
    },
}
with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
    json.dump(ora_config, tmp, indent=4)
    config_path = tmp.name

subprocess.run([binary, config_path], check=True)
os.remove(config_path)
print(f"  {dataset_name} -> {out_path}")