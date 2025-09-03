from goatools.utils import read_geneset
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
from goatools.base import download_ncbi_associations
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_go_basic_obo

# Load the GO DAG (ontology structure)
obo_basic = download_go_basic_obo()
obodag = GODag("go-basic.obo")

# annotations for genes from ncbi: gene2go
gene2go = download_ncbi_associations()
annotations = Gene2GoReader(gene2go, taxids=[9606], godag=obodag)
associations = annotations.get_ns2assc()

def calculate_terms(proximal_genes, background_genes):
    # Background genes: human coding genes
    background_genes = read_geneset(background_genes)

    gene_to_go = {}
    for ns, go_dict in associations.items():
        for gene, go_ids in go_dict.items():
            if gene not in gene_to_go:
                gene_to_go[gene] = set()
            gene_to_go[gene].update(go_ids)


    # GO enrichment study
    goea = GOEnrichmentStudy(
        background_genes,
        gene_to_go,
        obodag,
        methods=['fdr_bh']              # Multiple testing correction method (Benjamini-Hochberg FDR)
    )

    # GO enrichment analysis
    goea_results = goea.run_study(proximal_genes)

    significant_terms = []
    for r in goea_results:
        if r.p_fdr_bh < 0.05:
            significant_terms.append(r.GO)

    #goea.wr_xlsx("goea_face.xlsx", goea_results)

    return significant_terms


# Function to find all genes linked to a specific GO term
def find_genes_for_go_term(terms_list):
    geneid2gos = annotations.get_id2gos()  # Gene ID → GO Terms mapping

    go2gene = {}  # Dictionary to store GO term → Gene IDs
    for go_term_id in terms_list:
        go_term_id = go_term_id.strip()  # Remove extra characters/newlines

        if go_term_id not in obodag:
            print(f"GO term {go_term_id} not found in GO DAG")
            continue  # Skip missing GO terms

        # Get GO term and its descendants
        term = obodag[go_term_id]
        related_terms = term.get_all_children() | {go_term_id}

        # Find genes associated with these GO terms
        genes_linked = {gene for gene, gos in geneid2gos.items() if related_terms & gos}

        # Store results
        go2gene[go_term_id] = list(genes_linked)  # Convert set to list for JSON compatibility

    return go2gene



# calculate_terms(gene_IDs, "background_amiGO2_ids.txt")
