import csv
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Define your studies
studies = ["Breast_Cancer", "IBD", "T2D", "Schizophrenia", "Alzheimer"] #, "RA",
           #"Lupus", "Asthma", "MS", "Parkinson", "Prostate_Cancer"]

results = {
    "Study": [],
    "Category": [],
    "Count": []
}

for study in studies:
    folder = os.path.join("../eQTL", study)

    try:
        # Load gene sets
        with open(os.path.join(folder, "eqtl_genes.txt")) as f:
            eqtl_genes = set(line.strip() for line in f if line.strip())

        with open(os.path.join(folder, "GREAT_genes.txt")) as f:
            great_genes = set(line.strip() for line in f if line.strip())

        with open(os.path.join(folder, "novel_genes.txt")) as f:
            novel_genes = set(line.strip() for line in f if line.strip())

        # Compute counts
        results["Study"].extend([study] * 3)
        results["Category"].extend(["eQTL genes", "eQTL ∩ GREAT", "eQTL ∩ E-MAGO"])
        results["Count"].extend([
            len(eqtl_genes),
            len(eqtl_genes & great_genes),
            len(eqtl_genes & novel_genes)
        ])
        print(eqtl_genes & great_genes)
        print(eqtl_genes & novel_genes)


    except FileNotFoundError as e:
        print(f"Missing file for {study}: {e}")
        continue

# --- Plotting ---
df = pd.DataFrame(results)

plt.figure(figsize=(14, 6))
sns.barplot(data=df, x="Study", y="Count", hue="Category")
plt.xticks(rotation=45)
plt.title("Gene Set Overlaps Across 5 GWAS Studies")
plt.tight_layout()
plt.show()
