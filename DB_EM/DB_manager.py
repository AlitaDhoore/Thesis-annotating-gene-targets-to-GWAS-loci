# Class to manage the database and its connection
import mygene
import mysql.connector
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
import math

class DB_manager():
    def __init__(self):
        # connect
        self.db = mysql.connector.connect(
            host="localhost",
            user="root",
            password=""
        )

    def load_all_genes(self, background, found_genes):
        background_genes = set()
        if background is not None:
            with open(background, 'r') as file:
                for line in file:
                    background_genes.add(line.strip())
        all_genes = background_genes.union(found_genes)
        return all_genes

    def gene_exists(self, gene_name=None, gene_id=None):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        if gene_name:
            cursor.execute("SELECT gene_ID FROM gene_ID_name WHERE gene_name = %s", (gene_name,))
            result = cursor.fetchone()
            cursor.close()
            return result is not None
        elif gene_id:
            cursor.execute("SELECT gene_name FROM gene_ID_name WHERE gene_ID = %s", (gene_id,))
            result = cursor.fetchone()
            cursor.close()
            return result is not None
        cursor.close()
        return False

    # insert gene IDs and names into table and convert to each other
    # Function to convert gene ID to gene name
    def get_gene_name_from_id(self, gene_id):
        # Create a MyGeneInfo client
        mg = mygene.MyGeneInfo()
        result = mg.query(gene_id, scopes="entrezgene", fields="symbol", species="human")
        if "hits" not in result or not result["hits"]:
            print(f"No results found for gene {gene_id}")
            return None
        gene_name = result["hits"][0].get("symbol", "Not Found")
        if gene_name == "Not Found":
            print(f"for {gene_id} no gene name was found, gene is skipped.")
            return None
        return gene_name

    # Function to convert gene name to gene ID
    def get_gene_id_from_name(self, gene_name):
        # Create a MyGeneInfo client
        mg = mygene.MyGeneInfo()
        result = mg.query(gene_name, scopes="symbol", fields="entrezgene", species="human")
        if "hits" not in result or not result["hits"]:
            print(f"No results found for gene {gene_name}")
            return None
        gene_id = result["hits"][0].get("entrezgene", "Not Found")
        if gene_id == "Not Found":
            print(f"for {gene_name} no gene ID was found, gene is skipped.")
            return None
        return gene_id

    def insert_gene(self, gene):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        gene = str(gene)
        if not gene[0].isdigit():  # If it's a gene name
            gene_name = gene
            if self.gene_exists(gene_name=gene_name):  # Skip if already in DB
                return None
            gene_id = self.get_gene_id_from_name(gene_name)
        else:  # if it's a numeric gene IDD
            gene_id = int(gene)
            if self.gene_exists(gene_id=gene_id):  # Skip if already in DB
                return None
            gene_name = self.get_gene_name_from_id(gene_id)

        # Insert only if both values are present and update when a gene name or gene ID changes over time
        if gene_name and gene_id:
            cursor.execute("""
                INSERT INTO gene_ID_name (gene_name, gene_ID)
                VALUES (%s, %s)
                ON DUPLICATE KEY UPDATE gene_name = VALUES(gene_name), gene_ID = VALUES(gene_ID)
            """, (gene_name, gene_id))

        self.db.commit()
        cursor.close()

    def add_genes(self, background_file=None, found_genes=None):
        if found_genes is None:
            found_genes = []
        all_genes = self.load_all_genes(background_file, found_genes)

        for gene in all_genes:
            self.insert_gene(gene)

    def map_snps_to_tads(self, bed_file):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        # Ensure SNP_in_TAD table is empty before inserting data
        cursor.execute("SELECT COUNT(*) FROM SNP_in_TAD")
        count = cursor.fetchone()[0]

        if count != 0:
            cursor.execute("TRUNCATE TABLE SNP_in_TAD")

        # Read SNPs from the BED file
        with open(bed_file, 'r') as file:
            snps = [line.strip().split() for line in file]

        insert_query = """
            INSERT INTO SNP_in_TAD (SNP_ID, chromosome, Start_TAD, End_TAD)
            VALUES (%s, %s, %s, %s);
        """
        inside = 0
        between = 0
        outside = 0
        for chrom, _, pos, snp_id in snps:
            pos = int(pos)  # Convert position to integer

            # Fetch TADs for the chromosome
            cursor.execute("""
                SELECT Start_TAD, End_TAD 
                FROM TAD_Hg38
                WHERE chromosome = %s
                ORDER BY Start_TAD ASC;
            """, (chrom,))
            tads = cursor.fetchall()  # Fetch all results

            if not tads:
                print(f"Warning: No TADs found for chromosome {chrom}, SNP {snp_id}")
                continue

            previous_tad = None
            left_tad = None
            right_tad = None

            for start, end in tads:
                if start <= pos <= end:  # SNP inside a TAD
                    inside += 1
                    cursor.execute(insert_query, (snp_id, chrom, start, end))
                    break
                elif previous_tad and previous_tad[1] < pos < start:  # SNP between TADs
                    between += 1
                    left_tad = previous_tad[0]
                    right_tad = end
                    cursor.execute(insert_query, (snp_id, chrom, left_tad, right_tad))
                    break
                previous_tad = (start, end)

            # Handle SNP before the first TAD
            if pos < tads[0][0]:
                outside += 1
                cursor.execute(insert_query, (snp_id, chrom, tads[0][0], tads[0][1]))

            # Handle SNP after the last TAD
            elif pos > tads[-1][1]:
                outside += 1
                cursor.execute(insert_query, (snp_id, chrom, tads[-1][0], tads[-1][1]))
        print(f"inside = {inside}, between = {between}, outside = {outside}")
        self.db.commit()
        cursor.close()

    # this method links the genes to the SNPs according to the TAD they are situated in
    # SNPs without any genes linked --> save in another list --> take margin of 1 Mb
    def all_genes_linked_to_snps(self, margin=1000000):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        # Step 1: Clear SNP2gene table
        cursor.execute("SELECT COUNT(*) FROM SNP2gene")
        count = cursor.fetchone()[0]

        if count != 0:
            cursor.execute("TRUNCATE TABLE SNP2gene")

        # Step 2: Insert SNP-gene pairs linked by TAD (TAD_status = 1)
        cursor.execute("""
            INSERT INTO SNP2gene (SNP_ID, gene_ID, TAD_status)
            SELECT st.SNP_ID, gi.gene_ID, 1
            FROM SNP_in_TAD st
            JOIN coding_genes_in_TAD gt 
                ON st.chromosome = gt.chromosome 
                AND (
                    (st.Start_TAD = gt.Start_TAD AND st.End_TAD = gt.End_TAD) OR
                    (st.Start_TAD <= gt.Start_TAD AND st.End_TAD >= gt.End_TAD)
                )
            JOIN gene_ID_name gi 
                ON gt.Gene_Name = gi.Gene_Name;
                        """)

        # Step 3: Get SNPs that didn't get mapped to genes via TADs (TADless SNPs)
        cursor.execute("""
            SELECT DISTINCT st.SNP_ID, st.chromosome, CAST(SUBSTRING_INDEX(st.SNP_ID, ':', -1) AS UNSIGNED) AS snp_pos
            FROM SNP_in_TAD st
            LEFT JOIN SNP2gene sg
                ON st.SNP_ID = sg.SNP_ID
            WHERE sg.SNP_ID IS NULL
        """)

        tadless_snps = cursor.fetchall()

        print(f"Found {len(tadless_snps)} TADless SNPs")

        # Step 4: For each TADless SNP, find nearby genes in ±1Mb window
        insert_query = """
            INSERT IGNORE INTO SNP2gene (SNP_ID, gene_ID, TAD_status)
            VALUES (%s, %s, 0)
        """

        for snp_id, chrom, snp_pos in tadless_snps:
            cursor.execute("""
                SELECT gi.gene_ID
                FROM coding_genes_in_TAD cg
                JOIN gene_ID_name gi ON cg.gene_name = gi.gene_name
                WHERE cg.chromosome = %s
                  AND (
                      cg.Start_gene BETWEEN %s AND %s
                      OR cg.End_gene BETWEEN %s AND %s
                      OR (%s BETWEEN cg.Start_gene AND cg.End_gene)
                  );
            """, (chrom, snp_pos - margin, snp_pos + margin, snp_pos - margin, snp_pos + margin, snp_pos))

            nearby_genes = cursor.fetchall()

            for (gene_id,) in nearby_genes:
                cursor.execute(insert_query, (snp_id, gene_id))

        self.db.commit()

        cursor.execute("SELECT DISTINCT gene_ID FROM SNP2gene")
        all_genes = [row[0] for row in cursor.fetchall()]

        cursor.close()

        return all_genes

    def fetch_nearest_gene_to_SNP(self, bed_file):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        # Fetch all SNPs and their positions
        with open(bed_file, 'r') as file:
            snps = [line.strip().split() for line in file]

        # Dictionary to store closest gene per SNP
        closest_genes = {}

        for chrom, _, pos, snp_id in snps:
            pos = int(pos)  # Convert SNP position to integer

            # Fetch all gene IDs linked to this SNP from SNP2gene
            cursor.execute("""
                SELECT sg.gene_ID
                FROM SNP2gene sg
                WHERE sg.SNP_ID = %s;
            """, (snp_id,))

            gene_ids = [row[0] for row in cursor.fetchall()]

            if not gene_ids:
                print(f"No linked genes found for SNP: {snp_id}")
                continue

            # Retrieve gene positions from coding_genes_in_TAD by linking through gene_ID_name
            placeholders = ', '.join(['%s'] * len(gene_ids))  # Create placeholders for each gene_id
            query = f"""
                SELECT cg.gene_name, cg.Start_gene, cg.End_gene
                FROM coding_genes_in_TAD cg
                JOIN gene_ID_name gi ON cg.gene_name = gi.gene_name
                WHERE gi.gene_ID IN ({placeholders}) AND cg.chromosome = %s;
            """
            cursor.execute(query, (*gene_ids, chrom))

            rows = cursor.fetchall()

            if not rows:
                print(f"No gene positions found for SNP: {snp_id}")
                continue

            # Find the closest gene for the current SNP
            min_distance = float('inf')
            closest_gene_id = None

            for gene_name, start_gene, end_gene in rows:
                # Calculate distance to the gene
                distance = min(abs(pos - start_gene), abs(pos - end_gene))

                if distance < min_distance:
                    min_distance = distance

                    # Retrieve gene ID from gene_ID_name
                    cursor.execute("""
                        SELECT gene_ID FROM gene_ID_name WHERE gene_name = %s;
                    """, (gene_name,))

                    closest_gene_id = cursor.fetchone()[0]

            if closest_gene_id:
                closest_genes[snp_id] = closest_gene_id

        cursor.close()

        return list(closest_genes.values())

    def calculate_terms(self, candidate_genes, obodag, associations):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        # Background genes: human coding genes
        cursor.execute("SELECT DISTINCT gene_ID FROM gene_ID_name")
        background_genes = [row[0] for row in cursor.fetchall()]

        # Prepare gene to GO mapping
        gene_to_go = {}
        for ns, go_dict in associations.items():
            for gene, go_ids in go_dict.items():
                gene_to_go.setdefault(gene, set()).update(go_ids)

        # GO enrichment study
        goea = GOEnrichmentStudy(
            background_genes,
            gene_to_go,
            obodag,
            methods=['fdr_bh']
        )

        # Run enrichment analysis
        goea_results = goea.run_study(candidate_genes)

        # Sort by FDR-corrected p-value
        sorted_results = sorted(goea_results, key=lambda r: r.p_fdr_bh)

        significant_terms = []
        go_terms_info = []
        extra_terms_needed = 10

        # First collect all significant terms
        for r in sorted_results:
            if r.p_fdr_bh < 0.05:
                significant_terms.append(r.GO)
                go_terms_info.append({
                    "GO_ID": r.GO,
                    "FDR_p_value": r.p_fdr_bh,
                    "Description": r.name,
                    "studied_genes": r.study_items,
                    "Rank": len(go_terms_info) + 1
                })

        # If fewer than 10, add top non-significant ones
        if len(significant_terms) < 10:
            for r in sorted_results:
                if r.GO in significant_terms:
                    continue  # skip duplicates
                significant_terms.append(r.GO)
                go_terms_info.append({
                    "GO_ID": r.GO,
                    "FDR_p_value": r.p_fdr_bh,
                    "Description": r.name,
                    "studied_genes": r.study_items,
                    "Rank": len(go_terms_info) + 1
                })
                if len(significant_terms) >= 10:
                    break

            if len([r for r in go_terms_info if r['FDR_p_value'] < 0.05]) == 0:
                print("!Warning! No significant terms (FDR < 0.05) found. Top 10 terms used.")
            else:
                print("Note: Fewer than 10 significant terms found. Supplemented with top non-significant terms.")

        cursor.close()
        return significant_terms, go_terms_info

    def find_genes_for_go_terms(self, terms_list, obodag, annotations):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        new_terms_list = []
        for go_term_id in terms_list:
            cursor.execute("SELECT 1 FROM GO2gene WHERE GO_ID = %s LIMIT 1", (go_term_id,))
            exists = cursor.fetchone()
            if not exists:
                new_terms_list.append(go_term_id)
        if len(new_terms_list) != 0:
            geneid2gos = annotations.get_id2gos()  # Gene ID → GO Terms mapping

            for go_term_id in new_terms_list:
                go_term_id = go_term_id.strip()  # Remove extra characters/newlines

                if go_term_id not in obodag:
                    print(f"GO term {go_term_id} not found in GO DAG")
                    continue  # Skip missing GO terms

                # Get GO term and its descendants
                term = obodag[go_term_id]
                related_terms = term.get_all_children() | {go_term_id}

                # Find genes associated with these GO terms
                genes_linked = {gene for gene, gos in geneid2gos.items() if related_terms & gos}

                # Store results
                for linked_gene in genes_linked:
                    cursor.execute("""INSERT INTO GO2gene (GO_ID, gene_ID)
                    VALUES (%s, %s)
                    ON DUPLICATE KEY UPDATE GO_ID = GO_ID """, (go_term_id, linked_gene))

        self.db.commit()
        cursor.close()

    def fetch_genes_from_go_terms(self, terms_list):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        # Prepare a query to fetch gene IDs and their names related to GO terms
        format_strings = ','.join(['%s'] * len(terms_list))

        cursor.execute(f"""
            SELECT DISTINCT gene_ID_name.gene_name
            FROM GO2gene
            JOIN gene_ID_name ON GO2gene.gene_ID = gene_ID_name.gene_ID
            WHERE GO2gene.GO_ID IN ({format_strings})
        """, tuple(terms_list))

        # Fetch all gene names
        results = cursor.fetchall()

        # Extract gene names from the results and make them unique
        unique_genes = list(set(row[0] for row in results))

        cursor.close()
        return unique_genes

    def update_gene_scores(self, significant_terms, GO_terms_info):
        cursor = self.db.cursor(buffered=True)
        cursor.execute("USE SNP_DB")

        if significant_terms:
            format_strings = ','.join(['%s'] * len(significant_terms))

            # Convert the list of GO term information to a dictionary for easy lookup
            GO_term_scores = {
                go_info["GO_ID"]: 1 - math.log10(go_info["FDR_p_value"])
                for go_info in GO_terms_info if go_info["GO_ID"] in significant_terms
            }

            # Sum all significant GO term associations for each gene
            query = f"""
                UPDATE SNP2gene 
                JOIN (
                    SELECT gene_ID, SUM(COALESCE(gene_counts.count_score, 0)) AS total_score
                    FROM GO2gene
                    LEFT JOIN (
                        SELECT GO_ID, COUNT(*) AS count_score
                        FROM GO2gene
                        WHERE GO_ID IN ({format_strings})
                        GROUP BY GO_ID
                    ) AS gene_counts ON GO2gene.GO_ID = gene_counts.GO_ID
                    GROUP BY gene_ID
                ) AS gene_total 
                ON SNP2gene.gene_ID = gene_total.gene_ID
                SET SNP2gene.gene_score = COALESCE(SNP2gene.gene_score, 0) + gene_total.total_score
            """

            cursor.execute(query, tuple(significant_terms))

            # Update gene scores based on GO term significance
            for go_id, score in GO_term_scores.items():
                cursor.execute("""
                    UPDATE SNP2gene 
                    JOIN GO2gene ON SNP2gene.gene_ID = GO2gene.gene_ID
                    SET SNP2gene.gene_score = COALESCE(SNP2gene.gene_score, 0) + %s
                    WHERE GO2gene.GO_ID = %s
                """, (score, go_id))

        self.db.commit()
        cursor.close()

    def fetch_max_genes(self):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        cursor.execute("""
            SELECT SNP_ID, gene_ID, gene_score
            FROM (
                SELECT SNP_ID, gene_ID, gene_score, 
                       DENSE_RANK() OVER (PARTITION BY SNP_ID ORDER BY gene_score DESC) AS rnk
                FROM SNP2gene
            ) ranked
            WHERE rnk = 1;
        """)

        max_genes = cursor.fetchall()
        cursor.close()
        return max_genes

    def fetch_top_genes(self, bed_file):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")

        # Fetch top genes with a score higher than 0
        cursor.execute("""
            SELECT SNP_ID, gene_name, gene_score
            FROM (
                SELECT SNP2gene.SNP_ID, gene_ID_name.gene_name, SNP2gene.gene_score,
                       ROW_NUMBER() OVER (PARTITION BY SNP2gene.SNP_ID ORDER BY SNP2gene.gene_score DESC, gene_ID_name.gene_name ASC) AS rnk
                FROM SNP2gene
                JOIN gene_ID_name ON SNP2gene.gene_ID = gene_ID_name.gene_ID
                WHERE SNP2gene.gene_score > 0 
            ) ranked
            WHERE rnk <= 5;
        """)

        results = cursor.fetchall()

        # Convert results to dictionary {SNP_ID: [(gene_name, gene_score)]}
        snp_to_gene_data = {}
        for snp_id, gene_name, gene_score in results:
            if snp_id not in snp_to_gene_data:
                snp_to_gene_data[snp_id] = []
            snp_to_gene_data[snp_id].append((gene_name, gene_score))

        # Fetch all SNPs that do not have any gene annotations with score > 0
        cursor.execute("""
            SELECT DISTINCT sg.SNP_ID
            FROM SNP2gene sg
            LEFT JOIN gene_ID_name gi ON sg.gene_ID = gi.gene_ID
            WHERE sg.gene_score <= 0;
        """)
        zero_score_snps = [row[0] for row in cursor.fetchall()]

        # Fetch SNP positions from the BED file
        with open(bed_file, 'r') as file:
            snps = {line.strip().split()[3]: (line.strip().split()[0], int(line.strip().split()[2])) for line in file}

        for snp_id in zero_score_snps:
            if snp_id in snp_to_gene_data:
                continue  # Skip if already annotated with a valid gene

            if snp_id not in snps:
                continue  # Skip if SNP is not found in the BED file

            chrom, snp_pos = snps[snp_id]

            # Retrieve all genes linked to this SNP
            cursor.execute("""
                SELECT sg.gene_ID, cg.Start_gene, cg.End_gene, cg.gene_name
                FROM SNP2gene sg
                JOIN gene_ID_name gi ON sg.gene_ID = gi.gene_ID
                JOIN coding_genes_in_TAD cg ON gi.gene_name = cg.gene_name
                WHERE sg.SNP_ID = %s AND cg.chromosome = %s;
            """, (snp_id, chrom))

            rows = cursor.fetchall()

            if not rows:
                snp_to_gene_data[snp_id] = [("Intergenic", 0)]
                continue

            # Find the closest gene
            min_distance = float('inf')
            closest_gene = None

            for gene_id, start_gene, end_gene, gene_name in rows:
                distance = min(abs(snp_pos - start_gene), abs(snp_pos - end_gene))

                if distance < min_distance:
                    min_distance = distance
                    closest_gene = gene_name

            if closest_gene:
                snp_to_gene_data[snp_id] = [(closest_gene, 0)]

        # Fetch all SNPs that do not have any gene annotations at all
        cursor.execute("""
            SELECT DISTINCT st.SNP_ID
            FROM SNP_in_TAD st
            LEFT JOIN SNP2gene sg ON st.SNP_ID = sg.SNP_ID
            WHERE sg.SNP_ID IS NULL;
        """)
        geneless_snps = [row[0] for row in cursor.fetchall()]

        for snp_id in geneless_snps:
            if snp_id not in snp_to_gene_data:
                snp_to_gene_data[snp_id] = [("Intergenic", 0)]

        cursor.close()

        return snp_to_gene_data

    def monitor_gene_evolution(self):
        cursor = self.db.cursor()
        cursor.execute("USE SNP_DB")
        # Step 1: Query SNP2gene for rs10022462
        cursor.execute("SELECT gene_ID, gene_score FROM SNP2gene WHERE SNP_ID = 'rs71338792'")
        gene_rows = cursor.fetchall()

        gene2score = {}
        # Step 2: Loop through genes and extract gene_name + GO term count
        for gene_id, gene_score in gene_rows:
            # Get gene_name
            cursor.execute("SELECT gene_name FROM gene_ID_name WHERE gene_ID = %s", (gene_id,))
            result = cursor.fetchone()
            gene_name = result[0] if result else "Unknown"
            gene2score[gene_name] = gene_score, gene_id
        return gene2score







