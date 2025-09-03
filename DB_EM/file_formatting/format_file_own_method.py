import pandas as pd

# Load your TSV file
file_path = "../eQTL/ASD/ASD_GWAScatalog_leadSNPS.tsv"
df = pd.read_csv(file_path, sep="\t")

# Filter rows with valid 'locations' (contains a colon)
valid_locations = df[df['locations'].str.contains(":")]

# Split 'locations' into chromosome and position
location_split = valid_locations['locations'].str.split(':', expand=True)

# Clean positions: take the first number if there's a comma (just in case)
position_clean = location_split[1].str.split(',', expand=True)[0]

# Add 'chr' prefix to chromosome values
chromosomes = "chr" + location_split[0]

# Convert position to integers
positions = position_clean.astype(int)

# Calculate pos - 1
positions_minus_1 = positions - 1

# SNP IDs from the riskAllele column
snp_ids = valid_locations['riskAllele']

# Combine into a final dataframe
final_df = pd.DataFrame({
    'chromosome': chromosomes,
    'pos-1': positions_minus_1,
    'pos': positions,
    'SNP_ID': snp_ids
})

# Save without headers or index
output_path = "../eQTL/ASD/processsed_SNPs_ASD_GWAScatalog.txt"
final_df.to_csv(output_path, sep='\t', index=False, header=False)
