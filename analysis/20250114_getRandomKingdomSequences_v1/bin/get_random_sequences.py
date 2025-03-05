#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Francesco Costa 2025-02-07 fcosta@ebi.ac.uk

"""

import sys
sys.path.append("../../bin")
from pfamenv import PFAM_USER, PFAM_HOST, PFAM_PASSWORD, PFAM_PORT, PFAM_VERSION
from mysql import connector
import pandas as pd
from Bio import SeqIO

AFDB_SEQUENCES = "/nfs/research/agb/research/francesco/projects/20240212_isopeptideBonds_v1/20240529_findWithJess_v1/analysis/20240902_AFDBtoFasta_v1/output/AFDB.fa"
# Number of seqs to initilally extract from pfamseq
PFAM_NUMBER_SEQS=10000
# Number of sequences that should match to AFDB
AFDB_SEQS=1000

def main():
    
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)

    # Get random Bacteria sequences
    cursor = cnx.cursor()
    cursor.execute(f"SELECT pfamseq_acc,length,taxonomy\
                    FROM {PFAM_VERSION}.pfamseq\
                    WHERE taxonomy LIKE 'Bacteria;%'\
                    AND length <= 1280\
                    ORDER BY RAND()\
                    LIMIT {PFAM_NUMBER_SEQS};")
    output = cursor.fetchall()
    bac_df = pd.DataFrame(output, columns=["pfamseq_acc", "length", "taxonomy"])

    # Get random Archaea sequences
    cursor = cnx.cursor()
    cursor.execute(f"SELECT pfamseq_acc,length,taxonomy\
                    FROM {PFAM_VERSION}.pfamseq\
                    WHERE taxonomy LIKE 'Archaea;%'\
                    AND length <= 1280\
                    ORDER BY RAND()\
                    LIMIT {PFAM_NUMBER_SEQS};")
    output = cursor.fetchall()
    arc_df = pd.DataFrame(output, columns=["pfamseq_acc", "length", "taxonomy"])

    arc_count, bac_count = 0,0
    arc_ids, bac_ids = arc_df["pfamseq_acc"].unique(), bac_df["pfamseq_acc"].unique()

    for i in SeqIO.parse(AFDB_SEQUENCES, "fasta"):
        if i.id in arc_ids:
            arc_count += 1
        elif i.id in bac_ids:
            bac_count += 1
        if arc_count >= 1000 and bac_count >= 1000:
            break

    bac_df[bac_df["pfamseq_acc"].isin(bac_ids)].to_csv("output/bacteria_pfamseq_random.csv", index=False)
    arc_df[arc_df["pfamseq_acc"].isin(arc_ids)].to_csv("output/archaea_pfamseq_random.csv", index=False)

if __name__ == "__main__":
    main()