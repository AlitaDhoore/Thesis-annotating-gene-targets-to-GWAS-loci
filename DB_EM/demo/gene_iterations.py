import matplotlib.pyplot as plt

data_by_iteration = {}
current_iter = None

# === Read file ===
with open("snp2gene_iterations", "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("Iteration"):
            current_iter = int(line.split()[1])
            data_by_iteration[current_iter] = []
        elif line and current_iter is not None:
            parts = line.split("\t")
            if len(parts) >= 2:
                gene_name = parts[0]
                gene_score = int(parts[1])
                data_by_iteration[current_iter].append((gene_name, gene_score))

# === Generate horizontal barplots ===
for iteration, values in data_by_iteration.items():
    if not values:
        continue

    # Sort genes by score descending
    sorted_values = sorted(values, key=lambda x: x[1], reverse=True)
    genes, scores = zip(*sorted_values)

    plt.figure(figsize=(10, 6))
    plt.barh(genes, scores, color='skyblue')
    plt.xlabel("Gene Score")
    plt.title(f"Gene Scores - Iteration {iteration}")
    plt.gca().invert_yaxis()  # Highest score at the top
    plt.tight_layout()
    plt.show()

