import os
from goatools.utils import read_geneset
from goatools.obo_parser import GODag
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.go_enrichment import GOEnrichmentStudy

if "snakemake" not in dir():
    _ROOT = "/home/lukas/PycharmProjects/Ontologizer-Utils"
    _SET  = "GSE1297"
    _GEO  = f"{_ROOT}/resources/geo2kegg/{_SET}"

    class snakemake:
        config = {
            "go_obo": f"{_ROOT}/resources/go-basic.obo",
        }
        class input:
            study      = f"{_GEO}/{_SET}_study.txt"
            population = f"{_GEO}/{_SET}_population.txt"
            tab        = f"{_ROOT}/resources/goa_human_symbol.tab"
        output = [f"{_GEO}/{_SET}_freq_goatools_ora.csv"]

obo = snakemake.config["go_obo"]
tab = snakemake.input.tab

print("Loading ontology...")
godag = GODag(obo)

print("Loading annotations...")
id2gos = IdToGosReader(tab, godag=godag).get_id2gos()

out_path     = snakemake.output[0]
dataset_dir  = os.path.dirname(out_path)
dataset_name = os.path.basename(dataset_dir)

study_ids = read_geneset(snakemake.input.study)
pop_ids   = read_geneset(snakemake.input.population)

goeaobj = GOEnrichmentStudy(
    pop_ids, id2gos, godag,
    methods=["bonferroni", "fdr_bh"],
    pvalcalc="fisher_scipy_stats",
)
results = goeaobj.run_study_nts(study_ids)

with open(out_path, "w") as f:
    f.write("namespace,term_id,enrichment,p_uncorrected,p_fdr_bh,p_bonferroni,study_ratio,population_ratio\n")
    for ntd in sorted(results, key=lambda nt: [nt.p_uncorrected, nt.GO]):
        if ntd.p_fdr_bh < 0.05:
            f.write("{},{},{},{:.4e},{:.4e},{:.4e},{}/{},{}/{}\n".format(
                ntd.NS, ntd.GO, ntd.enrichment,
                ntd.p_uncorrected, ntd.p_fdr_bh, ntd.p_bonferroni,
                *ntd.ratio_in_study, *ntd.ratio_in_pop))

sig = sum(1 for ntd in results if ntd.p_fdr_bh < 0.05)
print(f"  {dataset_name}: {sig} significant terms -> {out_path}")