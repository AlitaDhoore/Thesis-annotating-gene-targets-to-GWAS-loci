import pandas as pd
# Replace these with your file paths
file1 = "IBD/gProfiler_hsapiens_29-05-2025_00-02-08__intersections.csv"
file2 = "IBD/gProfiler_novel_terms_IBD.csv"
file3 = "IBD/GREAT/gProfiler_terms_great.csv"
# great_file1 = "IBD/GREAT/GREAT_IBD_BP.tsv"
# great_file2 = "IBD/GREAT/GREAT_IBD_MF.tsv"
# great_file3 = "IBD/GREAT/GREAT_IBD_CC.tsv"
# Load the CSVs
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
df3 = pd.read_csv(file3)

# combine GREAT ids
# ids1 = pd.read_csv(great_file1, sep="\t", comment="#", header=None, usecols=[1, 5], names=["ID", "BinomBonfP"])
# ids2 = pd.read_csv(great_file2, sep="\t", comment="#", header=None, usecols=[1, 5], names=["ID", "BinomBonfP"])
# ids3 = pd.read_csv(great_file3, sep="\t", comment="#", header=None, usecols=[1, 5], names=["ID", "BinomBonfP"])


# Combine both GREAT outputs
# df_combined = pd.concat([ids1, ids2, ids3], ignore_index=True)

# Rank the dataframe based on "BinomBonfP" (low p-values to high p-values)
# df_combined_sorted = df_combined.sort_values(by='BinomBonfP', ascending=True)

# Ensure 'term_id' column exists
if "term_id" not in df1.columns or "term_id" not in df2.columns:
    raise ValueError("Both CSV files must contain a 'term_id' column.")


import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Define your sets
GREAT_terms = set(df3["term_id"])
novel_terms = set(df2["term_id"])
gc_terms = set(df1["term_id"])

gc_termids = df1["term_id"]
gc_termids.to_csv("IBD/terms_comparing/gc_terms.csv", index=False)


# Unique to Novel and GC (but not in GREAT)
novel_gc_unique = (novel_terms & gc_terms) - GREAT_terms

# Unique to GREAT and GC (but not in Novel)
great_gc_unique = (GREAT_terms & gc_terms) - novel_terms

# Unique to GC and Novel (but not in GREAT)
gc_great_novel = (gc_terms & GREAT_terms & novel_terms)

# Only keep 'source' and 'term_name' columns
unique_gc_novel = df1[df1['term_id'].isin(novel_gc_unique)][['term_id']]
unique_gc_great = df1[df1['term_id'].isin(great_gc_unique)][['term_id']]
common_all_three = df1[df1['term_id'].isin(gc_great_novel)][['term_id']]

unique_gc_novel.to_csv("IBD/terms_comparing/unique_gc_novel_terms.csv", index=False)
unique_gc_great.to_csv("IBD/terms_comparing/unique_gc_great_terms.csv", index=False)
common_all_three.to_csv("IBD/terms_comparing/common_all_three_terms.csv", index=False)


# Create the Venn diagram
plt.figure(figsize=(8, 8))
venn = venn3([gc_terms, novel_terms, GREAT_terms],
             set_labels=('GC', 'E-MAGO', 'GREAT'))

plt.title("Inflammatory Bowel Disease")
plt.show()