from goatools.anno.gaf_reader import GafReader

if "snakemake" not in dir():
    _ROOT = "/home/lukas/PycharmProjects/Ontologizer-Benchmark"

    class snakemake:
        class input:
            gaf = f"{_ROOT}/resources/goa_human.gaf"
        class output:
            tab = f"{_ROOT}/resources/goa_human_symbol.tab"

gaf_path = snakemake.input.gaf
tab_path = snakemake.output.tab

print(f"Building symbol tab from {gaf_path} -> {tab_path}")
rdr = GafReader(gaf_path)
gene2gos = {}
for rec in rdr.associations:
    gene = rec.DB_Symbol
    gene2gos.setdefault(gene, set()).add(rec.GO_ID)

with open(tab_path, "w") as f:
    for gene, gos in sorted(gene2gos.items()):
        f.write(f"{gene}\t{';'.join(sorted(gos))}\n")

print(f"  Wrote {len(gene2gos)} genes.")