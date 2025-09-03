import csv
import ast
import re
from EM_DB import EM_SNP_annotation

# Dictionary of genes with SNPs
snp2gene_great = {}
with open("../T2D/GREAT/snp2gene_GREAT_T2D.txt") as file:
    tsv_file = csv.reader(file, delimiter="\t")
    next(tsv_file, None)  # Skip header

    for line in tsv_file:
        if not line or len(line) < 2:
            continue

        snp = line[0].strip()
        gene_info = line[1].strip()
        if gene_info.startswith("NONE"):
            continue

        # Extract (gene, distance) pairs like HLA-DQA2 (-27636)
        gene_distance_pairs = []
        matches = re.findall(r'([\w\-]+)\s*\(([+-]?\d+)\)', gene_info)
        for gene, distance in matches:
            gene_distance_pairs.append((gene, abs(int(distance))))

        # Sort by absolute distance and store only gene names
        sorted_genes = tuple(g for g, _ in sorted(gene_distance_pairs, key=lambda x: x[1]))
        snp2gene_great[snp] = sorted_genes

GREAT_genes = []
for snp in snp2gene_great:
    GREAT_genes.append(snp2gene_great[snp][0])

result, info, memory = EM_SNP_annotation("T2D/processed_leadSNPs_T2D.txt", "../background_amiGO2_ids.txt")

import matplotlib.pyplot as plt

percentage_overlap = []

for i in range(len(memory)):
    same_genes = set(GREAT_genes).intersection(set(memory[i]))
    percent = (len(same_genes) / len(GREAT_genes)) * 100
    percentage_overlap.append(percent)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(range(len(percentage_overlap)), percentage_overlap, marker='o', linestyle='-', color='mediumblue')
plt.title('Type 2 Diabetes')
plt.xlabel('Number of Iterations')
plt.ylabel('Percentage of genes overlapping (%)')
plt.ylim(0, 100)
plt.grid(True)
plt.tight_layout()
plt.show()
