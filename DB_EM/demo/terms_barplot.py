import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim

# --- Load your file ---
file_path = "terms_iterations"

# Step 1: Parse terms per iteration
iterations = defaultdict(list)
current_iter = None

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("Iteration"):
            current_iter = line
        elif line.startswith("GO:") and current_iter:
            go_id = line.split()[0]
            iterations[current_iter].append(go_id)

go_dag = GODag("../go-basic.obo")
goslim_dag = GODag("../goslim_generic.obo")

# Step 2: Map and plot per iteration
for iter_name, go_ids in iterations.items():
    # Filter out invalid GO terms
    valid_terms = [go for go in go_ids if go in go_dag]

    # Map to slim
    all_slim = []
    for go in valid_terms:
        try:
            direct, _ = mapslim(go, go_dag, goslim_dag)
            all_slim.extend(direct)
        except ValueError:
            continue

    if not all_slim:
        print(f"{iter_name}: No GO Slim mappings found.")
        continue

    # Count frequencies
    slim_counts = Counter(all_slim)
    slim_labels = [goslim_dag[go].name if go in goslim_dag else go for go in slim_counts.keys()]
    counts = list(slim_counts.values())

    # Sort by count
    sorted_pairs = sorted(zip(slim_labels, counts), key=lambda x: x[1], reverse=True)
    labels, values = zip(*sorted_pairs)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.barh(labels, values, color='skyblue')
    plt.xlabel("Frequency")
    plt.title(f"GO Slim Term Distribution - {iter_name}")
    plt.gca().invert_yaxis()  # Highest count at top
    plt.tight_layout()
    plt.show()
