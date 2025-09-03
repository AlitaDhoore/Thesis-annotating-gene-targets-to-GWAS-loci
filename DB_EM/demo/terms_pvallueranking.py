import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict

# Load file
file_path = "terms_iterations"

# Step 1: Parse data per iteration
iterations = defaultdict(list)
current_iter = None

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("Iteration"):
            current_iter = line
        elif line.startswith("GO:") and current_iter:
            parts = line.split("\t")
            if len(parts) >= 3:
                go_id = parts[0]
                try:
                    pval = float(parts[1])
                except ValueError:
                    continue
                description = parts[2]
                iterations[current_iter].append((go_id, pval, description))

# Step 2: Plot top 10 per iteration by p-value
for iter_name, term_list in iterations.items():
    if not term_list:
        continue

    # Sort by p-value (lowest first)
    top_terms = sorted(term_list, key=lambda x: x[1])[:20]
    descriptions = [t[2] for t in top_terms]
    pvals = [t[1] for t in top_terms]
    neg_log_pvals = [-np.log10(p) if p > 0 else 0 for p in pvals]

    # Plot
    plt.figure(figsize=(10, 6))
    bars = plt.barh(descriptions, neg_log_pvals, color="salmon")
    plt.xlabel("-log10(FDR p-value)")
    plt.title(f"Top 20 GO Terms by FDR (lowest p-values) - {iter_name}")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()
