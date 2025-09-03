import mysql.connector
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

conn = mysql.connector.connect(
    host="localhost",
    user="root",
    password=""
)

cursor = conn.cursor()
cursor.execute("USE SNP_DB")

# # Query the data
# query = """
# SELECT SNP_ID, COUNT(DISTINCT gene_ID) AS num_genes
# FROM SNP2gene
# WHERE TAD_status = 1
# GROUP BY SNP_ID;
# """
# df1 = pd.read_sql(query, conn)

# Query the data
query = """
SELECT SNP_ID, COUNT(DISTINCT gene_ID) AS num_genes
FROM SNP2gene
WHERE TAD_status = 0
GROUP BY SNP_ID;
"""
df2 = pd.read_sql(query, conn)

# Close the connection
conn.close()
print(df2["num_genes"])
# Plot distribution
plt.figure(figsize=(10, 6))
sns.histplot(df2["num_genes"], bins=range(3, 15), kde=False)
plt.xlabel("Number of Genes per SNP")
plt.ylabel("Frequency")
plt.title("Distribution of Linked Genes per SNP +- 1 Mb")
plt.show()

