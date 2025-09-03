import mysql.connector
from mygene import MyGeneInfo

# Initialize MyGeneInfo
mg = MyGeneInfo()

# Connect to MySQL
db = mysql.connector.connect(
    host="localhost",
    user="root",
    password=""
)

cursor = db.cursor()

# Create database
cursor.execute("CREATE DATABASE IF NOT EXISTS SNP_DB")

# Use the database
cursor.execute("USE SNP_DB")

# Create tables
# Table to store TAD boundaries
cursor.execute("""
CREATE TABLE IF NOT EXISTS TAD_Hg38 (
    chromosome VARCHAR(10),
    Start_TAD INT,
    End_TAD INT,
    PRIMARY KEY (chromosome, Start_TAD, End_TAD)
);
""")

# Check if the table is empty
cursor.execute("SELECT COUNT(*) FROM TAD_Hg38")
count = cursor.fetchone()[0]  # Fetch the count result

if count == 0:  # If table is empty, insert data
    with open('../HiC_GM12878_DI_50kb.txt', 'r') as file:
        for line in file:
            columns = line.strip().split(' ')

            # Insert the data into the table
            cursor.execute("""
                INSERT INTO TAD_Hg38 (chromosome, Start_TAD, End_TAD)
                VALUES (%s, %s, %s)
            """, (columns[0], int(columns[1]), int(columns[2])))

# Table for fast conversions between different gene names and IDs
cursor.execute("""
CREATE TABLE IF NOT EXISTS gene_ID_name (
    gene_name VARCHAR(100),
    gene_ID INT,
    PRIMARY KEY (gene_name, gene_ID)
);
""")

# Table that has all coding genes per TAD
# useful to connect SNPs in TADs --> that way we know the coding genes in TADs
cursor.execute("""
CREATE TABLE IF NOT EXISTS coding_genes_in_TAD (
    gene_name VARCHAR(100),
    Start_gene INT,
    End_gene INT,
    chromosome VARCHAR(10),
    Start_TAD INT,
    End_TAD INT,
    PRIMARY KEY (gene_name, chromosome, Start_gene, End_gene)
);
""")


def fetch_genes_in_boundaries(chromosome, start, end, genome_version="hg38"):
    db_ucsc = mysql.connector.connect(
        host="genome-mysql.soe.ucsc.edu",
        user="genome",
        database=genome_version
    )

    # SQL query to fetch distinct genes within the specified boundaries
    # fetches genes that are fully in the TAD or partially
    query = f"""
    SELECT DISTINCT rg.name2 AS gene, rg.txStart, rg.txEnd
    FROM {genome_version}.refGene AS rg
    JOIN {genome_version}.wgEncodeGencodeAttrsV47 AS ga
    ON rg.name2 = ga.geneName
    WHERE rg.chrom = '{chromosome}' 
    AND (
        (rg.txStart >= {start} AND rg.txEnd <= {end})  
        OR (rg.txStart BETWEEN {start} AND {end})      
        OR (rg.txEnd BETWEEN {start} AND {end})      
    )
    AND ga.geneType = 'protein_coding'; 
    """

    cursor_ucsc = db_ucsc.cursor()
    cursor_ucsc.execute(query)
    genes = cursor_ucsc.fetchall()

    db_ucsc.close()

    # Remove duplicates by using a dictionary
    unique_genes = {gene[0]: {"found_gene": gene[0], "txStart": gene[1], "txEnd": gene[2]} for gene in genes}

    return list(unique_genes.values())


# Check if the table is empty
cursor.execute("SELECT COUNT(*) FROM coding_genes_in_TAD")
count = cursor.fetchone()[0]  # Fetch the count result
i = 0
if count == 0:
    cursor.execute("SELECT chromosome, Start_TAD, End_TAD FROM TAD_Hg38;")
    tad_boundaries = cursor.fetchall()

    for chromosome, start_tad, end_tad in tad_boundaries:
        genes_in_tad = fetch_genes_in_boundaries(chromosome, start_tad, end_tad)
        i += 1
        print(f"row: {i}, findings: {genes_in_tad}")
        for gene_data in genes_in_tad:
            cursor.execute("""
                INSERT INTO coding_genes_in_TAD (gene_name, Start_gene, End_gene, chromosome, Start_TAD, End_TAD)
                VALUES (%s, %s, %s, %s, %s, %s)
                ON DUPLICATE KEY UPDATE Start_gene = VALUES(Start_gene), 
                                        End_gene = VALUES(End_gene);
            """, (gene_data["found_gene"], gene_data["txStart"], gene_data["txEnd"], chromosome, start_tad, end_tad))

# Table to store input SNP into the right TAD
# if SNP not in TAD boundaries take the TADs before and behind it
cursor.execute("""
CREATE TABLE IF NOT EXISTS SNP_in_TAD (
    SNP_ID VARCHAR(50) PRIMARY KEY,
    chromosome VARCHAR(10),
    Start_TAD INT,
    End_TAD INT
);
""")

# Specify the coding genes present in the TADs where the SNPs reside
cursor.execute("""
CREATE TABLE IF NOT EXISTS SNP2gene (
    SNP_ID VARCHAR(50),
    gene_ID INT,
    gene_score INT DEFAULT 0,
    TAD_status TINYINT(1),
    PRIMARY KEY (SNP_ID, gene_ID)
);
""")

# save all used GO term to genes conversions
cursor.execute("""
CREATE TABLE IF NOT EXISTS GO2gene (
    GO_ID VARCHAR(50),
    gene_ID INT,
    PRIMARY KEY (GO_ID, gene_ID)
);
""")


def create_index_if_not_exists(cursor, table_name, index_name, column_name):
    # Check if the index already exists
    cursor.execute("""
        SELECT COUNT(1)
        FROM information_schema.statistics
        WHERE table_schema = DATABASE()
          AND table_name = %s
          AND index_name = %s;
    """, (table_name, index_name))

    exists = cursor.fetchone()[0]

    if not exists:
        # Create the index since it doesn't exist
        cursor.execute(f"CREATE INDEX {index_name} ON {table_name} ({column_name});")


create_index_if_not_exists(cursor, 'GO2gene', 'idx_GO2gene_GOID', 'GO_ID')
create_index_if_not_exists(cursor, 'SNP2gene', 'idx_SNP2gene_SNPID', 'SNP_ID')


db.commit()
cursor.close()
db.close()
