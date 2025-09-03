from gprofiler import GProfiler

# Initialize tests_alternatives
gp = GProfiler(return_dataframe=True)

with open('../genes.txt', 'r') as file:
    lines = file.readlines()

lines = [line.strip() for line in lines]

# Perform GO enrichment analysis
results = gp.profile(
    query=lines,
    organism="hsapiens",
    sources=['GO:BP', 'GO:MF', 'GO:CC'],
    no_evidences=False
)

results.to_excel('gprofiler_results.xlsx', index=False)
