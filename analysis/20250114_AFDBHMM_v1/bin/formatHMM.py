#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Format results table, add domain annotation
    Francesco Costa 2024-10-23 fcosta@ebi.ac.uk
    `sbatch --job-name="format_HMM" -t 48:00:00 --mem 16GB -e log/format_%j.err\
                 -o log/format_%j.out --mail-type=END\
                 --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/formatHMM.py"`

"""
import sys
sys.path.append("../../bin")
from pfamenv import PFAM_HOST, PFAM_PASSWORD, PFAM_PORT, PFAM_USER, PFAM_VERSION
from mysql import connector
import pandas as pd
import os
from Bio.SearchIO import HmmerIO

# Input
HMM = "output/AFDB.hmm"
# Output
OUTPUT_TABLE = "output/formatted_HMM.csv"
TAX_CACHE = "tmp/hmm_taxonomy.csv"
DOM_CACHE = "tmp/hmm_domains.csv"

def main():

    # Get results
    sequence_evalue = 0.01
    hit_evalue = 0.03
    outlist = []
    with open(HMM, "r") as fh:
        parser = HmmerIO.Hmmer3DomtabHmmqueryParser(fh)
        results = [line for line in parser.__iter__()]
        for hmm_ in results:
            for sequence in hmm_:
                if sequence.evalue >= sequence_evalue:
                    continue
                for hit in sequence:
                    if hit.evalue >= hit_evalue:
                        continue
                    outlist.append([hit.query_id, hit.hit_id, hit.query_range[0], hit.query_range[1],
                                   hit.hit_range[0], hit.hit_range[1], sequence.evalue, hit.evalue, hit.bitscore])
    hmm_df = pd.DataFrame(outlist, columns=["HMM", "uniprot_acc", "hmm_start", "hmm_end", "seq_start", "seq_end",
                                  "seq_evalue", "hit_evalue", "hit_bitscore"])
    hmm_df["Cna"] = hmm_df["HMM"].apply(lambda x: x.split("_")[5])

    # Map to taxonomy
    if os.path.isfile(TAX_CACHE) == False:
        tax_df = taxonomy(hmm_df["uniprot_acc"].unique())
        tax_df.to_csv(TAX_CACHE, index=False)
    else:
        tax_df = pd.read_csv(TAX_CACHE)
    
    tax_df["kingdom"] = tax_df["taxonomy"].apply(lambda x: x.split(";")[0].replace(".",""))
    
    # Map to domains
    if os.path.isfile(DOM_CACHE) == False:
        domains_df = domain(hmm_df["uniprot_acc"].unique())
        domains_df.to_csv(DOM_CACHE, index=False)
    else:
        domains_df = pd.read_csv(DOM_CACHE)
    domains_df = domains_df.rename(columns={"seq_start":"pfam_seq_start", "seq_end":"pfam_seq_end", "pfamseq_acc":"uniprot_acc"})

    # Add taxonomy
    hmm_df = pd.merge(hmm_df, tax_df, how="left").fillna("ND")
    
    # Add domains
    def isDomain(row):
        """
    
            Assign if at least 2/3 residues are in domain
        
        """
        seq_start = row["seq_start"]
        seq_end = row["seq_end"]
        pfam_seq_start = row["pfam_seq_start"]
        pfam_seq_end = row["pfam_seq_end"]
    
        matched_domain_len = seq_end-seq_start
        overlap = 0
        for aa in range(seq_start, seq_end):
            if aa >= pfam_seq_start and aa < pfam_seq_end:
                overlap += 1
        status = False
        
        # Overlap should cover at least 2/3 of the matched domain
        if overlap / matched_domain_len >= 2 / 3:
            status = True
        return status
    # Assign domains and delete unnecessary rows 
    hmm_df = pd.merge(hmm_df, domains_df[["uniprot_acc", "pfamA_acc", "pfam_seq_start", "pfam_seq_end", "pfamA_id"]], how="left")
    hmm_df["is_domain"] = hmm_df.apply(lambda x: isDomain(x), axis=1)
    hmm_df = hmm_df.sort_values("is_domain", ascending=False)\
        .drop_duplicates(["HMM", "uniprot_acc", "seq_start", "seq_end"], keep="first")

    hmm_df.to_csv(OUTPUT_TABLE, index=False)


def domain(uniprot_accs:list) -> pd.DataFrame:
    """
    
        Maps uniprot accs to pfam domains
    
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    uniprot_accs = ",".join([f"'{u}'" for u in uniprot_accs])

    query = f"""SELECT uniprot_acc, uniprot_reg_full.pfamA_acc, pfamA.pfamA_id, clan_membership.clan_acc, clan.clan_id, seq_start, seq_end \
                FROM {PFAM_VERSION}.uniprot_reg_full, {PFAM_VERSION}.clan_membership, {PFAM_VERSION}.clan, {PFAM_VERSION}.pfamA
                WHERE uniprot_acc IN ({uniprot_accs}) \
                AND uniprot_reg_full.pfamA_acc=clan_membership.pfamA_acc \
                AND clan.clan_acc=clan_membership.clan_acc \
                AND uniprot_reg_full.pfamA_acc=pfamA.pfamA_acc\
                AND uniprot_reg_full.in_full=1"""

    cursor.execute(query)
    output = cursor.fetchall()

    return(pd.DataFrame(output, columns=["uniprot_acc", "pfamA_acc", "pfamA_id", "clan_acc", "clan_id", "seq_start", "seq_end"]))

def taxonomy(uniprot_accs:list) -> pd.DataFrame:
    """
    
         Maps uniprot accs to taxonomy
    
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    uniprot_accs = ",".join([f"'{u}'" for u in uniprot_accs])
    # For this use old pfam (closer to AFDB)
    cursor.execute(f"SELECT uniprot_acc, taxonomy, species FROM {PFAM_VERSION}.uniprot \
                WHERE uniprot_acc IN ({uniprot_accs})")

    output = cursor.fetchall()

    return(pd.DataFrame(output, columns=["uniprot_acc", "taxonomy", "species"]))

if __name__ == "__main__":
    main()