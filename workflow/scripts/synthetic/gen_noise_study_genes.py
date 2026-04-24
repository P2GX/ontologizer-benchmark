import os
from common import term_gene_map, generate_bench_data, write_bench_data

if "snakemake" not in dir():
    class snakemake:
        class input:
            go = "../../resources/go/go-basic.obo"
            gaf = "../../resources/go/goa_human.gaf"
        class params:
            recall = "0.4"
        class wildcards:
            noise = "0.2"
        output = ["../../resources/synthetic/noise/noise_0.2.tsv"]
    os.makedirs("../../resources/synthetic/noise", exist_ok=True)

go_path = snakemake.input.go
gaf_path = snakemake.input.gaf
out_path = snakemake.output[0]

noise = float(snakemake.wildcards.noise)
fixed_recall = float(snakemake.params.recall)

term_to_genes = term_gene_map(go_path, gaf_path)

print(f"Noise sweep: recall={fixed_recall}, noise={noise}")
repeat_data = {i: generate_bench_data(term_to_genes, recall=fixed_recall, noise=noise, n_genes=500)
               for i in range(10)}
write_bench_data(repeat_data, out_path)

