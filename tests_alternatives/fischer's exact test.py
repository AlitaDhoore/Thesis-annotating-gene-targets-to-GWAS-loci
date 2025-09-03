# we test if the genes annotated with the GO-term are overrepresented in the input set compared to the background set
# so for example if we have 98 genes in the input set and 60 of them are associated with GO term and 38 not
# we have 20,000 background genes where also 60 are associated with term

from scipy.stats import fisher_exact
res = fisher_exact([[60, 38], [60, 19840]], alternative='greater')
print(res.pvalue)

# we conclude that the odds of having a gene in the input set associated with the GO-term are high compared to
# the odds of a background gene being associated with the GO-term
# low pvalue means enrichment in term