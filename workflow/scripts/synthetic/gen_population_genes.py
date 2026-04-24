from common import term_gene_map

if "snakemake" not in dir():
    class snakemake:
        class input:
            go = "../../resources/go/go-basic.obo"
            gaf = "../../resources/go/goa_human.gaf"
        output = ["../../resources/synthetic/population_genes.txt"]

go_path = snakemake.input.go
gaf_path = snakemake.input.gaf
out_path = snakemake.output[0]

term_to_genes = term_gene_map(go_path, gaf_path)

population = set.union(*term_to_genes.values())
with open(out_path, 'w') as f:
    for gene in population:
        f.write(f"{gene}\n")