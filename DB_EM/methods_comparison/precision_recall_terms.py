import pandas as pd
import matplotlib.pyplot as plt

# Function to load term sets from a file
def load_terms(filepath):
    return set(pd.read_csv(filepath)["term_id"])

# Function to compute precision and recall
def compute_pr(predicted, truth):
    tp = len(predicted & truth)
    fp = len(predicted - truth)
    fn = len(truth - predicted)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    return precision, recall

# File paths for each disease
datasets = {
    "T2D": {
        "GC": "T2D/gProfiler_GC_terms_T2D.csv",
        "E-MAGO": "T2D/gProfiler_novel_terms_T2D.csv",
        "GREAT": "T2D/GREAT/gProfiler_terms_great.csv"
    },
    "IBD": {
        "GC": "IBD/gProfiler_GC_terms_IBD.csv",
        "E-MAGO": "IBD/gProfiler_novel_terms_IBD.csv",
        "GREAT": "IBD/GREAT/gProfiler_terms_great.csv"
    },
    "BC": {
        "GC": "breast_carcinoma/gProfiler_GO_terms_GC.csv",
        "E-MAGO": "breast_carcinoma/gProfiler_GO_terms_novel.csv",
        "GREAT": "breast_carcinoma/GREAT/gProfiler_terms_great.csv"
    }
}

# Colors and markers for plotting
colors = {
    "T2D": "royalblue",
    "IBD": "seagreen",
    "BC": "tomato"
}
markers = {
    "E-MAGO": "s",
    "GREAT": "o"
}

# Plotting
plt.figure(figsize=(8, 6))

for disease, files in datasets.items():
    gc_terms = load_terms(files["GC"])
    emago_terms = load_terms(files["E-MAGO"])
    great_terms = load_terms(files["GREAT"])

    prec_emago, rec_emago = compute_pr(emago_terms, gc_terms)
    prec_great, rec_great = compute_pr(great_terms, gc_terms)

    plt.scatter(rec_emago, prec_emago, color=colors[disease], marker=markers["E-MAGO"],
                label=f"{disease} - E-MAGO", s=100, edgecolor='black')
    plt.scatter(rec_great, prec_great, color=colors[disease], marker=markers["GREAT"],
                label=f"{disease} - GREAT", s=100, edgecolor='black')

plt.xlabel("Recall (vs GC)")
plt.ylabel("Precision (vs GC)")
plt.title("Precision-Recall for GO terms")
plt.xlim(0, 1.05)
plt.ylim(0, 1.05)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
