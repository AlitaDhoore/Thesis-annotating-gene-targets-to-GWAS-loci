import re

# --- Step 1: Parse GREAT file and collect distant SNPs
great_file = "Alzheimer/snp2gene_GREAT_alz.txt"
filtered_snps = set()

with open(great_file, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        snp = parts[0]
        genes = parts[1].split(", ")

        far = True
        for gene in genes:
            match = re.search(r"\(([+-]?)(\d+)\)", gene)
            if match:
                distance = int(match.group(2))
                if distance <= 50000:
                    far = False
                    break
        if far:
            filtered_snps.add(snp)

print(f"Retained {len(filtered_snps)} SNPs >50kb from nearest gene.")

# --- Step 2: Filter the SNP-based files with custom column index
def filter_snp_file(input_file, output_file, snp_set, snp_col_index=0):
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) > snp_col_index and parts[snp_col_index] in snp_set:
                    fout.write(line)


# File paths
file0 = "Alzheimer/processed_leadSNPs_alz.txt"   # SNP in column 3 (index 2)
file1 = "Alzheimer/snp_gene_alz.txt"             # SNP in column 1 (index 0)
file2 = "Alzheimer/snp2gene_GREAT_alz.txt" # SNP in column 1 (index 0)

# Apply filtering with appropriate SNP column
filter_snp_file(file0, "Alzheimer/lonely_SNPs.txt", filtered_snps, snp_col_index=3)
filter_snp_file(file1, "Alzheimer/lonely_novel.txt", filtered_snps, snp_col_index=0)
filter_snp_file(file2, "Alzheimer/lonely_great.txt", filtered_snps, snp_col_index=0)
