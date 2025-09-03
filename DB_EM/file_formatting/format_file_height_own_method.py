import pandas as pd

# Load the TSV file
file_path = "../eQTL/Alzheimer/leadSNPs_alz.tsv"
df = pd.read_csv(file_path, sep="\t")

# Ensure necessary columns exist and drop rows with missing values
df = df.dropna(subset=['CHR_ID', 'CHR_POS', 'SNPS'])

# Convert 'CHR_ID' to string and add 'chr' prefix (handles X, Y, MT, etc.)
df['CHR_ID'] = 'chr' + df['CHR_ID'].astype(float).astype(int).astype(str)

# Convert 'CHR_POS' to integer
df['CHR_POS'] = df['CHR_POS'].astype(int)

# Create 'pos-1' column (1-based position minus 1 for 0-based format)
df['POS_MINUS_1'] = df['CHR_POS'] - 1

# Select relevant columns and rename them
final_df = df[['CHR_ID', 'POS_MINUS_1', 'CHR_POS', 'SNPS']]
final_df.columns = ['chromosome', 'pos-1', 'pos', 'SNP_ID']

# Define output file path
output_path = "../eQTL/Alzheimer/processed_leadSNPs_alz.txt"

# Save the processed data without headers or index
final_df.to_csv(output_path, sep='\t', index=False, header=False)

print(f"Processed file saved as: {output_path}")

