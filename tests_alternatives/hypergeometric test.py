import numpy as np
from scipy.stats import hypergeom
import matplotlib.pyplot as plt

# every ontology annotation separately
# k successes = overlap of found genes with genes with ontology annotation
# K possible successes = genes with ontology annotation pi
# n number of draws = the genes found (after first iteration this is proximity)
# N population size = all genes with an ontology annotation

# variation in k
N = 10000
K = 500
n = 50

# Range of possible successes in the sample: 0 to min(K, n)
k = np.arange(0, min(K, n) + 1)

# Calculate the PMF for each value of x
result = hypergeom.pmf(k, N, K, n)

# Plot the PMF
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with N=10000')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()

# if k goes down the pvalue will decrease

N = 1000
result = hypergeom.pmf(k, N, K, n)
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with N=1000')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()

N = 100000
result = hypergeom.pmf(k, N, K, n)
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with N=100000')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()

# if N goes up the pvalue will increase

N = 1000
n = 100
result = hypergeom.pmf(k, N, K, n)
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with n=100')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()

N = 1000
n = 10
result = hypergeom.pmf(k, N, K, n)
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with n=10')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()

N = 1000
n = 50
K = 900
result = hypergeom.pmf(k, N, K, n)
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with K=900')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()

N = 1000
K = 100
result = hypergeom.pmf(k, N, K, n)
plt.bar(k, result, color='skyblue', edgecolor='black')
plt.title('Hypergeometric PMF with K=100')
plt.xlabel('Number of successes (k)')
plt.ylabel('Probability')
plt.show()