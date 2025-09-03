# Step 1: initialisation:
# input = BED file with lead SNPs + BED file with TAD boundaries
# Start with lead SNPs and find their TADs
    # TADKB uses Hg19 as reference genome --> if loop to see if lead SNPs use hg19 or hg 38 --> use liftover (UCSC) to convert

import pandas as pd
import mysql.connector
from mygene import MyGeneInfo
from liftover import get_lifter

mg = MyGeneInfo()


# find the TAD boundaries of the SNPs by looking through a text file from TADKB
def get_tad_boundaries(chromosome, position, boundaries):
    TADs_data = open(boundaries, 'r')
    previous_TAD = [0, 0, 0]

    TADs = TADs_data.readlines()
    for TAD in TADs:
        TAD = TAD.split()
        chrom = TAD[0]
        begin = TAD[1]
        end = TAD[2]
        if chrom == chromosome:
            if int(begin) <= int(position) <= int(end):
                return begin, end
            elif int(previous_TAD[1]) <= int(position) <= int(end):
                return previous_TAD[1], end
        previous_TAD = TAD


    return None


# converting with liftover tool --> we need to first check which genome it is!!!
def process_snps(input_bed, boundaries, TADs_file):
    snp_data = pd.read_csv(input_bed, sep='\t', header=None, names=["chrom", "pos-1", "pos", "snp_id"])

    boundaries_data = []

    # We need to convert the hg38 to hg19 because TADKB only offers boundaries for hg19
    # Loop through each SNP in the input data
    for index, row in snp_data.iterrows():
        chrom = row["chrom"]
        pos = row["pos"]
        snp_id = row["snp_id"]

        # we need to convert hg19 to hg38 with liftover tool
        # converter = get_lifter('hg38', 'hg19', one_based=True)
        # new_pos = converter[chrom][pos]
        # pos = new_pos[0][1]

        # Use the start position for querying the TAD boundary
        # if boundaries can be found append them to boundary list
        if get_tad_boundaries(chrom, pos, boundaries) is None:
            print(f"Skipping {snp_id} due to missing TAD boundaries.")
        else:
            tad_start, tad_end = get_tad_boundaries(chrom, pos, boundaries)
            boundaries_data.append([chrom, tad_start, tad_end, snp_id])

        # we need to convert hg19 to hg38 with liftover tool
        # converter = get_lifter('hg19', 'hg38', one_based=True)
        # new_tad_start = converter[chrom][int(tad_start)]
        # new_tad_end = converter[chrom][int(tad_end)]
        # tad_start = new_tad_start[0][1]
        # tad_end = new_tad_end[0][1]


    # Write to output BED file
    output_df = pd.DataFrame(boundaries_data, columns=["chrom", "tad_start", "tad_end", "snp_id"])
    output_df.to_csv(TADs_file, sep='\t', header=False, index=False)


# Within these TADs we need to summarize all genes present using UCSC MySQL database
def fetch_genes_in_boundaries(chromosome, start, end, genome_version="hg38"):
    db = mysql.connector.connect(
        host="genome-mysql.soe.ucsc.edu",
        user="genome",
        database=genome_version
    )

    # SQL query to fetch genes within the specified boundaries
    query = f"""
    SELECT rg.name2 AS gene, rg.txStart, rg.txEnd
    FROM {genome_version}.refGene AS rg
    JOIN {genome_version}.wgEncodeGencodeAttrsV47 AS ga
    ON rg.name2 = ga.geneName
    WHERE rg.chrom = '{chromosome}' 
    AND rg.txStart >= {start} 
    AND rg.txEnd <= {end}
    AND ga.geneType = 'protein_coding'; 
    """

    cursor = db.cursor()
    cursor.execute(query)
    genes = cursor.fetchall()

    db.close()

    # Extract gene names, start of gene and end of gene from the query result
    return [{"found_gene": gene[0], "txStart": gene[1], "txEnd": gene[2]} for gene in genes]


def find_genes_per_SNP(input_bed, boundaries, TADs_file, genes_file, genome_version="hg38"):
    process_snps(input_bed, boundaries, TADs_file)
    boundaries_found = pd.read_csv(TADs_file, sep='\t', header=None, names=["chrom", "start", "end", "snp_id"])

    snp2gene = {}
    all_genes_found = []
    initial_genes = []

    # Loop through each TAD boundary and fetch genes
    for index, row in boundaries_found.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        snp_id = row["snp_id"]

        # Fetch genes within the specified boundary
        genes = fetch_genes_in_boundaries(chrom, start, end, genome_version)
        if len(genes) == 0:
            print(f"no genes found in this TAD from {start} to {end} containing snp {snp_id}")
        # Order the genes from closest to furthest from the snp
        genes_w_dist_to_snp = {}
        ordered_genes = []
        for gene in genes:
            snp_coord = int(snp_id.split(':')[1])
            dist_start = abs(snp_coord-gene["txStart"])
            dist_end = abs(snp_coord-gene["txEnd"])
            closest_end = min(dist_start, dist_end)
            genes_w_dist_to_snp[gene["found_gene"]] = closest_end
            genes_w_dist_to_snp = dict(sorted(genes_w_dist_to_snp.items(), key=lambda item: item[1]))
        for i in range(0, len(genes_w_dist_to_snp)):
            ordered_genes.append(list(genes_w_dist_to_snp)[i])

        gene_ids = []
        if len(ordered_genes) > 0:
            for i in range(0, len(ordered_genes)):
                result = mg.query(ordered_genes[i], scopes="symbol", fields="entrezgene", species="human")
                gene_id = result["hits"][0].get("entrezgene", "Not Found")
                if gene_id == "Not Found":
                    print(f"for {ordered_genes[i]} no gene ID was found, gene is skipped.")
                else:
                    gene_ids.append(gene_id)
                    initial_genes.append(int(gene_id))

        snp2gene[snp_id] = gene_ids

        # Join genes into a comma-separated string
        gene_list = ",".join(ordered_genes)

        # Append SNP and gene list to output data
        all_genes_found.append([chrom, start, end, snp_id, gene_list])

    # Write to output BED file with each SNP and its list of genes
    output_df = pd.DataFrame(all_genes_found, columns=["chrom", "start", "end", "snp_id", "genes"])
    output_df.to_csv(genes_file, sep='\t', header=False, index=False)
    return snp2gene, initial_genes
