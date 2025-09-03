import requests
import pandas as pd
import numpy as np

# Functions
def query_variant(var_id, p_upper=5e-2):
    print(f"Querying for {var_id}...")
    url = f"https://www.ebi.ac.uk/eqtl/api/associations/{var_id}"
    params = {
        "paginate": False,
        "p_upper": p_upper
    }
    response = requests.get(url, params=params)
    if response.status_code != 200:
        print(f"Request failed: {response.status_code}")
        return False, pd.DataFrame()

    data = response.json()

    # Check if '_embedded' is in the response
    if '_embedded' not in data or 'associations' not in data['_embedded']:
        print(f"No associations found for {variant_id}")
        return False, pd.DataFrame()

    data = response.json()['_embedded']['associations']
    return response.ok, pd.DataFrame.from_dict(data).T

# Read rsID
input_name = f"eQTL/Alzheimer/lonely_SNPs.txt"
snp_column = pd.read_table(input_name, sep="\t", usecols=[3], header=None)

variant_list = np.unique(snp_column[3])

# Loop over vars
full_df = pd.DataFrame()
for variant_id in variant_list:
    ok, df = query_variant(variant_id, p_upper=5e-2)
    if ok and not df.empty:
        relevant_terms = [
            "cortex",
            "hippocampus",
            "precuneus",
            "brain"
        ]
        pattern = "|".join(relevant_terms)

        if 'qtl_group' in df.columns:
            filtered_df = df[df['qtl_group'].str.contains(pattern, case=False, na=False)]
            if not filtered_df.empty:
                full_df = pd.concat([
                    full_df,
                    filtered_df[['rsid', 'gene_id', 'pvalue', 'qtl_group', 'study_id']]
                ], ignore_index=True)
        else:
            print(f"'qtl_group' column missing for {variant_id}")


# Save
output_name = f"eQTL/Alzheimer/new_eQTLs.txt"
full_df.to_csv(output_name,index=False)