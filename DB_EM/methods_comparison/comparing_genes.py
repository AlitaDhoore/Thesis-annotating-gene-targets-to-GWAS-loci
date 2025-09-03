import csv
import ast
import re

first_novel_genes = []
novel_genes = []
# with open("eQTL/T2D/novel_genes.txt") as file:
#     tsv_file = csv.reader(file, delimiter="\t")
#     for line in tsv_file:
#         first_novel_genes.append(line[1])
#         for i in range(1, len(line), 2):
#             novel_genes.append(line[i])
#
# with open("eQTL/T2D/novel_genes.txt", "w") as f:
#     for gene in first_novel_genes:
#         f.write(gene + "\n")

# dictionary of genes with snps
snp2gene_great = {}
with open("../IBD/GREAT/snp2gene_GREAT_IBD.txt") as file:
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

with open("../IBD/novel_genes.txt") as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for line in tsv_file:
        novel_genes.append(line[0])

GREAT_genes = []
for snp in snp2gene_great:
    GREAT_genes.append(snp2gene_great[snp])

GC_genes = []
with open("../IBD/GC_genes_IBD.txt") as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for line in tsv_file:
        GC_genes.append(line[0])

print(set(GC_genes).intersection(set(novel_genes)))
print(set(GC_genes).intersection(set(GREAT_genes)))

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Define your sets
gc = set(GC_genes)
novel = set(novel_genes)
great = set(GREAT_genes)

# Create the Venn diagram
plt.figure(figsize=(8, 8))
venn = venn3([gc, novel, great],
             set_labels=('GC', 'E-MAGO', 'GREAT'))

plt.title("Inflammatory Bowel Disease")
plt.show()
