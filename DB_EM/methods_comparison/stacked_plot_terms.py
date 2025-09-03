import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load all three files
df_a = pd.read_csv("../IBD/terms_comparing/go_slim_frequency_gc.csv")  # Replace with your file paths
df_b = pd.read_csv("../IBD/terms_comparing/go_slim_frequency_gcvsnovel.csv")
df_c = pd.read_csv("../IBD/terms_comparing/go_slim_frequency_gcvsgreat.csv")
df_d = pd.read_csv("../IBD/terms_comparing/go_slim_frequency_allthree.csv")


# Standardize column names if necessary
df_a.columns = ['GO_Slim_ID', 'GO_Slim_Name', 'Frequency']
df_b.columns = ['GO_Slim_ID', 'GO_Slim_Name', 'Frequency']
df_c.columns = ['GO_Slim_ID', 'GO_Slim_Name', 'Frequency']
df_d.columns = ['GO_Slim_ID', 'GO_Slim_Name', 'Frequency']

# Merge on GO_Slim_Name with outer join to keep all categories
merged = pd.merge(df_a, df_b, on='GO_Slim_Name', how='outer', suffixes=('_a', '_b'))
merged = pd.merge(merged, df_c, on='GO_Slim_Name', how='outer')
merged.rename(columns={'Frequency': 'Frequency_c'}, inplace=True)
merged = pd.merge(merged, df_d, on='GO_Slim_Name', how='outer')
merged.rename(columns={'Frequency': 'Frequency_d'}, inplace=True)

# Fill NaNs with 0 and convert to int
merged = merged.fillna(0)
freq_cols = ['Frequency_a', 'Frequency_b', 'Frequency_c', 'Frequency_d']
merged[freq_cols] = merged[freq_cols].astype(int)

# Calculate total frequency across all sets
merged['Total'] = merged[freq_cols].sum(axis=1)

# Filter out low-frequency categories
merged = merged[merged['Total'] >= 3]

# Sort by total frequency for clearer plotting
merged = merged.sort_values(by='Total', ascending=False)

# Plotting
labels = merged['GO_Slim_Name']
a_counts = merged['Frequency_a'].values
b_counts = merged['Frequency_b'].values
c_counts = merged['Frequency_c'].values
d_counts = merged['Frequency_d'].values

x = np.arange(len(labels))
bar_width = 0.6

plt.figure(figsize=(12, 7))
plt.bar(x, a_counts, bar_width, label='GC terms', color='skyblue')
plt.bar(x, b_counts, bar_width, bottom=a_counts, label='Overlapping terms GC and novel', color='salmon')
plt.bar(x, c_counts, bar_width, bottom=a_counts + b_counts, label='Overlapping terms GC and GREAT', color='lightgreen')
plt.bar(x, d_counts, bar_width, bottom=a_counts + b_counts + c_counts, label='Overlap of all three', color='plum')

plt.xticks(x, labels, rotation=90)
plt.ylabel('Frequency')
plt.title('Inflammatory Bowel Disease')
plt.legend()
plt.tight_layout()
plt.show()
