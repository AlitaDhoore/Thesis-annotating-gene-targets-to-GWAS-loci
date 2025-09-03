import pandas as pd
import csv
import ast
import re

# Load your CSV file
df = pd.read_csv("Alzheimer/new_eQTLs.txt")  # Replace with your actual filename

# Convert scientific notation to floats (if needed)
df["pvalue"] = pd.to_numeric(df["pvalue"], errors='coerce')

# Drop rows with missing p-values
df = df.dropna(subset=["pvalue"])

# Group by rsid and get the row with the minimum p-value per group
min_pval_genes = df.loc[df.groupby("rsid")["pvalue"].idxmin()]

# Extract the gene_id column
gene_list = min_pval_genes[["gene_id"]]

# Save to a new text file with one gene_id per line
gene_list.to_csv("eQTL/Alzheimer/eqtl_genes.txt", index=False, header=False)

first_novel_genes = []
novel_genes = []
with open("Alzheimer/lonely_novel.txt") as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for line in tsv_file:
        first_novel_genes.append(line[1])
        for i in range(1, len(line), 2):
            novel_genes.append(line[i])

with open("Alzheimer/novel_genes.txt", "w") as f:
    for gene in first_novel_genes:
        f.write(gene + "\n")

# dictionary of genes with snps
snp2gene_great = {}
with open("Alzheimer/lonely_great.txt") as file:
    tsv_file = csv.reader(file, delimiter="\t")
    next(tsv_file, None)  # Skip header

    for line in tsv_file:
        if len(line) == 0:
            continue

        snp = line[0].strip()
        gene_info = line[1].strip()
        if gene_info[0] == "NONE":
            continue

        # Extract (gene, distance) tuples
        gene_distance_pairs = []
        for g in gene_info.split(','):
            match = re.match(r'(\w+)\s+\(([+-]?\d+)\)', g.strip())
            if match:
                gene = match.group(1)
                distance = int(match.group(2))
                gene_distance_pairs.append((gene, abs(distance)))

        # Sort by absolute distance and store only gene names
        sorted_genes = tuple(g for g, _ in sorted(gene_distance_pairs, key=lambda x: x[1]))
        if len(sorted_genes) == 0:
            continue
        snp2gene_great[snp] = sorted_genes

GREAT_genes = []
with open("Alzheimer/GREAT_genes.txt", "w") as f:
    for snp in snp2gene_great:
        GREAT_genes.append(snp2gene_great[snp][0])
        f.write(snp2gene_great[snp][0] + "\n")

