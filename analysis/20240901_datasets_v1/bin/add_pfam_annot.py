#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Add Pfam classification
    Francesco Costa 2025-10-02 fcosta@ebi.ac.uk

"""


import sys
sys.path.append("../../bin")
from pfamenv import PFAM_HOST,PFAM_PASSWORD,PFAM_PORT,PFAM_USER,PFAM_VERSION
from mysql import connector
import numpy as np
import pandas as pd
import os
from dotenv import load_dotenv

load_dotenv("../../.env")
TABLE = os.getenv("TABLE")
OUTPUT = "output/table_with_pfam_annot.csv"

def main():

    df = pd.read_csv(TABLE)
    # Get pfam info associated with each PDB
    columns = ["pdb_id", "pfamA_acc", "pfamA_id", "chain", "pdb_res_start", "pdb_res_end", "seq_start", "seq_end", "clan_acc", "clan_id"]

    proteins = df["PDB code"].unique()

    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()

    proteins_ = ", ".join([f"'{protein}'" for protein in proteins])

    # Get domains
    cursor.execute(f"SELECT pdb_id,pfamA.pfamA_acc,pfamA.pfamA_id,chain,pdb_res_start,pdb_res_end,seq_start,seq_end,\
                        clan_membership.clan_acc,clan.clan_id \
                    FROM {PFAM_VERSION}.pdb_pfamA_reg, {PFAM_VERSION}.pfamA, \
                        {PFAM_VERSION}.clan_membership, {PFAM_VERSION}.clan \
                    WHERE pdb_id IN ({proteins_}) \
                    AND pfamA.pfamA_acc = clan_membership.pfamA_acc \
                    AND clan.clan_acc=clan_membership.clan_acc \
                    AND pfamA.pfamA_acc = pdb_pfamA_reg.pfamA_acc")
    output = cursor.fetchall()
    d_df = pd.DataFrame(output, columns = columns)
    d_df["pdb_id"] = d_df["pdb_id"].apply(lambda x: x.lower())
    data_df = pd.merge(df, d_df.rename(columns={"pdb_id":"PDB code", "chain":"Chain"}), how="left")\
            .rename(columns={"Position 1\r\n(Bond 1)":"r1_bond", "Position 2\r\n(catalytic)":"r_cat", "Position 3\r\n(Bond 2)":"r2_bond"})

    # Check which domains map to isopep bonds
    data_df["is_domain"] = data_df.apply(lambda x: is_domain(x), axis=1)
    data_df = data_df.sort_values("is_domain", ascending=False) \
            .drop_duplicates(["PDB code", "r1_bond", "r_cat", "r2_bond"], keep="first")
    data_df.loc[data_df["is_domain"]==False, "pfamA_acc"] = np.NaN
    data_df.loc[data_df["is_domain"]==False, "pfamA_id"] = np.NaN
    data_df.loc[data_df["is_domain"]==False, "clan_acc"] = np.NaN
    data_df.loc[data_df["is_domain"]==False, "clan_id"] = np.NaN
    data_df.loc[data_df["is_domain"]==False, "seq_start"] = np.NaN
    data_df.loc[data_df["is_domain"]==False, "seq_end"] = np.NaN
    data_df["PDB code"] = data_df["PDB code"].apply(lambda x: x.upper())
    data_df.to_csv(OUTPUT, index=False)

def is_domain(row):
    """

        Assign if at least 2/3 residues are in domain
    
    """
    residues = [row["r1_bond"], row["r_cat"], row["r2_bond"]]
    status = False
    c = 0
    for res in residues:
        if res >= row["seq_start"] and res <= row["seq_end"]:
            c += 1
            #status = True
            #break
    if c >= 2:
        status = True
    return status


if __name__ == "__main__":
    main()