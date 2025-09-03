# 1. initialisation
# initialise the genes annotated to the SNP --> done by proximity
# how many genes? how will the genes be represented?
# the parameters that need to be updated are the GO terms found with GOEA

# 2. Expectation
# "Derive the expectation of the complete log-likelihood"
# Involves finding the expected value of hidden data
# in this case: we have GO terms with p-values and there might be genes relevant to these GO terms but not yet present
# in the current gene set
# select genes: associations in gene ontology database --> take the genes significantly associated to these GO terms

# 3. Maximization
# "optimize parameters by differentiating the expected log-likelihood"
# in this case: with new gene set update GO terms with by rerunning GOEA

# 4. Repeat
# new GO terms are used to repeat the process until convergence

import GOEA
import initialisation

def EM_SNP_annotation(input_snps, background_genes):
    gene_scores = {}
    old_candidate_genes = []

    # the initialisation needs a bed file with lead SNPs and a bed file with TAD boundaries
    # the output is a list with all genes closest to the lead SNPs
    snp2genes, initial_genes = initialisation.find_genes_per_SNP(input_snps, "TADs_hg38.bed", "../SNPs_withTADs.bed", "genes_in_TADs.bed")
    new_candidate_genes = initial_genes

    # keeping track of the iterations
    i = 1
    # keeping track of convergence
    c = 0
    while c < 0.75*len(new_candidate_genes):
        print(f"the comparison of genes after iteration {i} are \n new: {new_candidate_genes} \n old: {old_candidate_genes}")
        i += 1

        # 1. do GOEA 2. fetch significant GO-terms
        significant_terms = GOEA.calculate_terms(new_candidate_genes, background_genes)
        print(f"the significant terms after GOEA {significant_terms} ")

        # 3. find all genes related to these GO terms
        terms2genes = GOEA.find_genes_for_go_term(significant_terms)

        # 4. Count the amount of times a certain gene is found by all the significant terms and give gene scores
        for term in terms2genes:
            for gene in terms2genes[term]:
                if gene not in gene_scores:
                    gene_scores[gene] = 1
                else:
                    gene_scores[gene] += 1

        # 5. per locus/SNP take the gene with the highest gene score
        old_candidate_genes = new_candidate_genes
        new_candidate_genes = []
        for snp in snp2genes:
            highest_score = 0
            for gene in snp2genes[snp]:
                if int(gene) in gene_scores:
                    if gene_scores[int(gene)] > highest_score:
                        highest_score = gene_scores[int(gene)]
            if highest_score == 0:
                for gene in snp2genes[snp]:
                    new_candidate_genes.append(gene)
            else:
                highest_scoring_gene = list(gene_scores.keys())[list(gene_scores.values()).index(highest_score)]
                new_candidate_genes.append(highest_scoring_gene)

        c = 0
        for j in range(0, len(new_candidate_genes)):
            if old_candidate_genes[j] == new_candidate_genes[j]:
                c += 1
        print(f"the current similarity is {c}")

        # now do GOEA again until convergence!
    return None


EM_SNP_annotation("peaks_face_v1_hg38.bed", "background_amiGO2_ids.txt")