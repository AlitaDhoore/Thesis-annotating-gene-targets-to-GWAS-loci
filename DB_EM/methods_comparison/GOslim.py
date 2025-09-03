from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim
import csv
import os

# Step 1: Load full GO DAG
print("Downloading and loading GO DAG...")
obo_file = download_go_basic_obo()
go_dag = GODag(obo_file)  # This must be a GODag instance
print(f"Loaded {len(go_dag)} GO terms.")

# Step 2: Extract GO Slim term IDs from OBO file (text-based parsing)
goslim_file = "../goslim_generic.obo"
if not os.path.exists(goslim_file):
    raise FileNotFoundError(f"GO Slim file '{goslim_file}' not found.")

goslim_dag = GODag(goslim_file)

# Step 3: Read GO terms from CSV
input_csv = "IBD/terms_comparing/gc_terms.csv"
go_terms = set()

print(f"Reading GO terms from {input_csv}...")
with open(input_csv, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)  # skip header
    for row in reader:
        if row and row[0].startswith("GO:"):
            go_terms.add(row[0].strip())

print(f"Read {len(go_terms)} unique GO terms.")

# Step 4: Warn about missing terms
missing_terms = [term for term in go_terms if term not in go_dag]
if missing_terms:
    print(f"Warning: These GO terms are not in the GO DAG and will be skipped:\n{missing_terms}")

# Filter out missing ones
go_terms = [term for term in go_terms if term in go_dag]
print(go_terms)

# Step 6: Perform the mapping
print("Mapping GO terms to GO Slim...")
mapped = {}
for go_term in go_terms:
    try:
        direct, all_ = mapslim(go_term, go_dag, goslim_dag)
        mapped[go_term] = direct
    except ValueError:
        print(f"Skipping GO term not in GO DAG: {go_term}")
# Step 7: Display results
print("\nGO Term Mappings:")
for go_id, slim_ids in mapped.items():
    print(f"{go_id} maps to: {', '.join(slim_ids) if slim_ids else 'No mapping'}")

# # Optional: Save to CSV
# output_csv = "mapped_to_slim.csv"
# print(f"\nSaving mapping to '{output_csv}'...")
# with open(output_csv, "w", newline='') as outcsv:
#     writer = csv.writer(outcsv)
#     writer.writerow(["Original_GO", "Mapped_Slim_Terms"])
#     for go_id, slim_ids in mapped.items():
#         writer.writerow([go_id, ", ".join(slim_ids) if slim_ids else "No mapping"])

from collections import Counter

# Step 1: Flatten all GO slim mappings into one list
all_mapped_slim_terms = [slim for slims in mapped.values() for slim in slims]

# Step 2: Count frequency of each GO Slim term
slim_term_counts = Counter(all_mapped_slim_terms)

# Step 3: Map GO Slim term ID to its name/description
slim_term_names = {go_id: goslim_dag[go_id].name for go_id in slim_term_counts if go_id in goslim_dag}

# Step 4: Combine ID, name, and count into one list
slim_summary = []
for go_id, count in slim_term_counts.items():
    name = slim_term_names.get(go_id, "Unknown")
    slim_summary.append((go_id, name, count))

# Step 5: Sort by frequency
slim_summary = sorted(slim_summary, key=lambda x: x[2], reverse=True)

# Optional: Print top 10
print("\nTop mapped GO Slim terms:")
for go_id, name, count in slim_summary[:10]:
    print(f"{go_id} - {name}: {count} times")

# Optional: Save to CSV
output_freq_csv = "IBD/terms_comparing/go_slim_frequency_gc.csv"
import csv
with open(output_freq_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["GO_Slim_ID", "GO_Slim_Name", "Frequency"])
    writer.writerows(slim_summary)

print(f"\nFrequency summary saved to {output_freq_csv}")
