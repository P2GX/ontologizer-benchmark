import os
from common import term_gene_map, generate_bench_data, write_bench_data

go_path = snakemake.input.go
gaf_path = snakemake.input.gaf
out_path = snakemake.output[0]

recall = float(snakemake.wildcards.recall)
fixed_precision = float(snakemake.params.precision)

term_to_genes = term_gene_map(go_path, gaf_path)

print(f"Recall sweep: recall={recall}, precision={fixed_precision}")
repeat_data = {i: generate_bench_data(term_to_genes, recall=recall, precision=fixed_precision, n_genes=500)
               for i in range(10)}
write_bench_data(repeat_data, out_path)


