import os
from common import term_gene_map, generate_bench_data, write_bench_data

if "snakemake" not in dir():
    class snakemake:
        class input:
            go = "../../resources/go/go-basic.obo"
            gaf = "../../resources/go/goa_human.gaf"
        class params:
            noise = "0.4"
        class wildcards:
            recall = "0.2"
        output = ["../../resources/synthetic/recall/recall_0.2.tsv"]
    os.makedirs("../../resources/synthetic/recall", exist_ok=True)

go_path = snakemake.input.go
gaf_path = snakemake.input.gaf
out_path = snakemake.output[0]

recall = float(snakemake.wildcards.recall)
fixed_noise = float(snakemake.params.noise)

term_to_genes = term_gene_map(go_path, gaf_path)

print(f"Recall sweep: recall={recall}, noise={fixed_noise}")
repeat_data = {i: generate_bench_data(term_to_genes, recall=recall, noise=fixed_noise, n_genes=500)
               for i in range(10)}
write_bench_data(repeat_data, out_path)


