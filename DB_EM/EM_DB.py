from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_go_basic_obo
from DB_manager import DB_manager
from goatools.base import download_ncbi_associations

import csv


def EM_SNP_annotation(input_snps, background_genes):
    # Preparing the variables for GOEA
    # Load the GO DAG (ontology structure)
    obo_basic = download_go_basic_obo()
    obodag = GODag(obo_basic)

    # annotations for genes from ncbi: gene2go
    gene2go = download_ncbi_associations()
    annotations = Gene2GoReader(gene2go, taxids=[9606], godag=obodag)
    associations = annotations.get_ns2assc()

    # first we add the background genes into the database
    # using the lead SNPs we put all the genes that are found to be in the same TAD as the SNP in the database
    db_connect = DB_manager()
    db_connect.add_genes(background_file=background_genes)
    db_connect.map_snps_to_tads(input_snps)

    used_genes = db_connect.all_genes_linked_to_snps()
    initial_genes = db_connect.fetch_nearest_gene_to_SNP(input_snps)

    db_connect.add_genes(found_genes=used_genes)

    new_candidate_genes = initial_genes
    old_candidate_genes = ["initialisation"]

    # keeping track of the iterations
    i = 1
    # keeping track of convergence
    c = 0
    GO_terms_info = []
    genes_memory = []
    candidate_gene_names = []
    # for gene in new_candidate_genes:
    #     candidate_gene_names.append(db_connect.get_gene_name_from_id(gene))
    # genes_memory.append(candidate_gene_names)
    while c < min(len(old_candidate_genes), len(new_candidate_genes)):
        # 1. do GOEA 2. fetch significant GO-terms
        significant_terms_results = db_connect.calculate_terms(new_candidate_genes, obodag, associations)
        significant_terms = significant_terms_results[0]
        GO_terms_info = significant_terms_results[1]
        print(f"the significant terms after GOEA {significant_terms} ")

        # 3. find all genes related to these GO terms
        db_connect.find_genes_for_go_terms(significant_terms, obodag, annotations)

        # 4. Count the amount of times a certain gene is found by all the significant terms and give gene scores
        # with DB loop over term, query the genes for this term and add score in SNP2gene table
        db_connect.update_gene_scores(significant_terms, GO_terms_info)

        # 5. per locus/SNP take the gene with the highest gene score
        old_candidate_genes = new_candidate_genes
        new_candidate_genes = []

        max_genes = db_connect.fetch_max_genes()
        print(f"max score genes: {max_genes} and participating snps is {len(max_genes)}")

        # Extract gene_IDs into a list
        new_candidate_genes = [row[1] for row in max_genes]

        c = 0
        for j in range(0, min(len(old_candidate_genes), len(new_candidate_genes))):
            if old_candidate_genes[j] == new_candidate_genes[j]:
                c += 1
        print(f"the current similarity is {c}")
        print(f"the comparison of genes after iteration {i} are \n new: {new_candidate_genes} \n old: {old_candidate_genes}")
        candidate_gene_names = []
        # for gene in new_candidate_genes:
        #     candidate_gene_names.append(db_connect.get_gene_name_from_id(gene))
        # genes_memory.append(candidate_gene_names)

        # with open("demo/terms_iterations", "a") as f:
        #     f.write("Iteration" + " " + str(i) + "\n")
        #     for info in GO_terms_info:
        #         f.write(str(info["GO_ID"]) + "\t" + str(info["FDR_p_value"]) + "\t" + str(info["Description"]) + "\n")
        # with open("demo/snp2gene_iterations", "a") as f:
        #     f.write("Iteration" + " " + str(i) + "\n")
        #     gene2score = db_connect.monitor_gene_evolution()
        #     for gene in gene2score:
        #         f.write(f"{gene}\t{gene2score[gene][0]}\n")
        i += 1

        # now do GOEA again until convergence!

    result = db_connect.fetch_top_genes(input_snps)

    # Create a list to hold each line of the file
    lines = []

    # Iterate through the SNP dictionary
    for snp_id, genes in result.items():
        # Flatten gene names and scores into alternating values
        gene_info = []
        for gene_name, gene_score in genes:
            gene_info.append(gene_name)
            gene_info.append(str(gene_score))  # convert score to string for writing

        # Combine SNP ID and gene info into a single line (tab-separated)
        line = "\t".join([snp_id] + gene_info)
        lines.append(line)


    # Write the lines to a text file
    output_file = "eQTL/Alzheimer/snp_gene_alz_modified.txt"
    with open(output_file, "w") as f:
        for line in lines:
            f.write(line + "\n")

    return result, GO_terms_info #, genes_memory


(EM_SNP_annotation("eQTL/Alzheimer/processed_leadSNPs_alz.txt","../background_amiGO2_ids.txt"))