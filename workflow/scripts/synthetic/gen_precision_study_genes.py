import os
from common import term_gene_map, generate_bench_data, write_bench_data

go_path = snakemake.input.go
gaf_path = snakemake.input.gaf
out_path = snakemake.output[0]

precision = float(snakemake.wildcards.precision)
fixed_recall = float(snakemake.params.recall)

term_to_genes = term_gene_map(go_path, gaf_path)

print(f"Precision sweep: recall={fixed_recall}, precision={precision}")
repeat_data = {i: generate_bench_data(term_to_genes, recall=fixed_recall, precision=precision, n_genes=500)
               for i in range(10)}
write_bench_data(repeat_data, out_path)