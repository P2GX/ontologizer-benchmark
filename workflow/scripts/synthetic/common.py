import random
from typing import Dict, Set
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader

def term_gene_map(go_path : str, gaf_path : str) -> Dict[str, Set[str]]:
    # Parses the OBO file to resolve hierarchy and alternate IDs
    go_graph: GODag = GODag(go_path)

    # Parses the GAF file, handling format variations and evidence codes
    gaf_reader: GafReader = GafReader(gaf_path)

    term_to_genes: Dict[str, Set[str]] = {}

    for association in gaf_reader.associations:
        gene_id: str = association.DB_Symbol
        term_id: str = association.GO_ID

        # Verify the term exists in the current ontology version
        if term_id not in go_graph:
            continue

        # Exclude NOT-qualified annotations (gene experimentally shown NOT to have this function)
        if association.Qualifier and 'NOT' in association.Qualifier:
            continue

        # Restrict to protein-coding genes
        if association.DB_Type != 'protein':
            continue

        ancestors = go_graph[term_id].get_all_parents()
        term_to_genes.setdefault(term_id, set()).add(gene_id)
        for ancestor in ancestors:
            term_to_genes.setdefault(ancestor, set()).add(gene_id)

    return term_to_genes


def recall_set(
        genes: Set[str],
        recall: float
) -> Set[str]:
    tp_genes = random.sample(list(genes), int(len(genes)*recall))
    return set(tp_genes)

def generate_bench_data(term_to_genes, recall: float, noise: float, n_genes: int) -> Dict[str, Set[str]]:
    all_terms = set(term_to_genes.keys())
    all_genes = set.union(*term_to_genes.values())

    eligible_terms = [t for t in all_terms if 10 <= len(term_to_genes[t]) <= 200]

    # Accumulate signal genes from true terms until n_genes is reached
    signal_term_to_genes: Dict[str, Set[str]] = {}
    n_signal = 0
    while n_signal < n_genes:
        term = random.choice(eligible_terms)
        recall_genes = recall_set(term_to_genes[term], recall)
        signal_term_to_genes[term] = recall_genes
        n_signal += len(recall_genes)

    # Add noise genes on top: noise_level * n_signal random genes from population
    result = dict(signal_term_to_genes)
    if noise > 0:
        result["Noise"] = set(random.sample(sorted(all_genes), int(n_signal * noise)))

    return result


def write_bench_data(repeat_data: Dict[int, Dict[str, Set[str]]], path: str) -> None:
    with open(path, 'w') as f:
        f.write("repeat\tterm_id\tgene\n")
        for repeat, term_to_genes in repeat_data.items():
            for term, genes in term_to_genes.items():
                for gene in genes:
                    f.write(f"{repeat}\t{term}\t{gene}\n")


def extract_study_genes(tsv_path: str, repeat: int, out_path: str) -> None:
    with open(tsv_path) as f:
        next(f)  # skip header
        genes = {line.split('\t')[2].strip() for line in f if int(line.split('\t')[0]) == repeat}
    with open(out_path, 'w') as f:
        for gene in genes:
            f.write(f"{gene}\n")