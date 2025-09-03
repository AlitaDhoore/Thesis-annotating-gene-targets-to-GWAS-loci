import matplotlib.pyplot as plt

# Function to load gene list from a text file
def load_genes(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip().upper() for line in f if line.strip())

# === MANUALLY specify file paths ===
datasets = [
    {
        'disease': 'BC',
        'gc_file': 'breast_carcinoma/genecards_genes.txt',
        'great_file': 'breast_carcinoma/GREAT/GREAT_genes.txt',
        'novel_file': 'breast_carcinoma/novel_genes_BC.txt'
    },
    {
        'disease': 'IBD',
        'gc_file': 'IBD/GC_genes_IBD.txt',
        'great_file': 'IBD/GREAT/GREAT_genes.txt',
        'novel_file': 'IBD/novel_genes.txt'
    },
    {
        'disease': 'T2D',
        'gc_file': 'T2D/GC_genes_T2D',
        'great_file': 'T2D/GREAT/GREAT_genes.txt',
        'novel_file': 'T2D/novel_genes.txt'
    }
]

# Define markers and colors
methods = {
    'GREAT': 'o',   # circle
    'E-MAGO': 's'    # square
}
colors = {
    'BC': 'tomato',
    'IBD': 'seagreen',
    'T2D': 'royalblue'
}

# Initialize plot
plt.figure(figsize=(8, 6))

# Loop over each dataset entry
for data in datasets:
    disease = data['disease']
    gc_genes = load_genes(data['gc_file'])

    for method, marker in methods.items():
        file_path = data['great_file'] if method == 'GREAT' else data['novel_file']
        pred_genes = load_genes(file_path)

        # Compute metrics
        true_pos = len(pred_genes & gc_genes)
        precision = true_pos / len(pred_genes) if pred_genes else 0
        recall = true_pos / len(gc_genes) if gc_genes else 0

        # Plot
        plt.scatter(
            recall,
            precision,
            label=f'{disease} - {method}',
            color=colors[disease],
            marker=marker,
            s=100,
            edgecolor='black'
        )

# Final plot adjustments
plt.xlabel('Recall (vs GC)', fontsize=12)
plt.ylabel('Precision (vs GC)', fontsize=12)
plt.title('Precision–Recall Plot for Genes', fontsize=14)
plt.grid(True)
plt.xlim(0, 0.2)
plt.ylim(0, 0.2)
plt.legend(title='Disease - Method', fontsize=10, title_fontsize=11, loc='upper right')
plt.tight_layout()
plt.show()
