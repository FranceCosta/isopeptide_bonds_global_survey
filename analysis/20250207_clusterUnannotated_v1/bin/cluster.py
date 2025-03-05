#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    1. Extract sequences which are unmapped in Pfam 
    with isopep bond prob > 0.65 and shorter than 
    300 in interbond residues distance (insertions)
    2. Cluster them
    3. Reformat the output table
    `conda activate mmseqs2`
    Francesco Costa 2025-02-07 fcosta@ebi.ac.uk

"""

import pandas as pd
import subprocess
import os
from dotenv import load_dotenv
load_dotenv("../../.env")
AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_SCAN_SEQUENCES = os.getenv("AFDB_SCAN_SEQUENCES")
AFDB_DOMAINS = os.getenv("AFDB_DOMAINS")
TABLE_WITH_PFAM = os.getenv("TABLE_WITH_PFAM")

def main():
    df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
    d_df = pd.read_csv(AFDB_DOMAINS)
    seq_df = pd.read_csv(AFDB_SCAN_SEQUENCES)
    pdf = pd.read_csv(TABLE_WITH_PFAM)

    # Simplify everything.
    # This is a very rought method but will help
    # Get isopep bonds with no annotations:
    cond1 = (df["probability"]>.65)
    cond2 = (~df["taxonomy"].isna())
    cond3 = (df["pfamA_acc"].isna())
    final_df = pd.merge(df[cond1&cond2&cond3], seq_df)

    # Create domain boundaries
    # Apply different boundaries based on average distance between bond and domain boundary
    boundaries = {
        "CnaB-like": {
            "left": 10,
            "right": 5
        },
        "CnaA-like": {
            "left": 5,
            "right": 30
        }
    }
    final_df["seq_start"] = final_df.apply(lambda x: x["r1_bond"] - boundaries[x["bond_type"]]["left"], axis=1)
    final_df.loc[final_df["seq_start"]< 0, "seq_start"] = 1
    final_df["seq_end"] = final_df.apply(lambda x: x["r2_bond"] + boundaries[x["bond_type"]]["right"], axis=1)
    final_df["id"] = final_df.apply(lambda x: f"{x['uniprot_acc']}_{x['seq_start']}_{x['seq_end']}", axis=1)
    final_df["sequence"] = final_df.apply(lambda x: x["sequence"][x["seq_start"]:x["seq_end"]], axis=1)

    # Exclude long proteins (insertions most likely)
    final_df["len"] = final_df["seq_end"] - final_df["seq_start"]
    final_df = final_df[final_df["len"]<300]

    # Save sequences
    with open("output/sequences.fa", "wt") as fh:
        for index, row in final_df.iterrows():
            fh.write(">{}\n{}\n".format(row["id"], row["sequence"]))

    # Run mmseqs2
    run_mmseqs_clustering("output/sequences.fa", "output/idp", "tmp")
    
    # Reformat output
    clus_df = pd.read_table("output/idp_cluster.tsv", names=["clus_rep", "members"])
    clus_df = pd.merge(clus_df, final_df[["id", "template"]].rename(columns={"id":"members"}))
    # Add pfam annot
    pdf["template"] = pdf.apply(lambda x: x["PDB code"].lower()+"_"+x["Chain"]+"_"+str(x["r1_bond"])+"_"+str(x["r_cat"])+"_"+str(x["r2_bond"]), axis=1)
    clus_df = pd.merge(clus_df, pdf[["template", "pfamA_acc"]], how="left", on="template")
    # Collapse members and pfam annots
    clus_df = pd.merge(  
        pd.merge(clus_df.groupby("clus_rep")["members"].apply(lambda x: ','.join(x)).reset_index(),
            clus_df.fillna("ND").rename(columns={"pfamA_acc":"members_template_domains"}).groupby("clus_rep")["members_template_domains"].apply(lambda x: ','.join(x)).reset_index()),
        clus_df[["clus_rep", "template", "pfamA_acc"]].rename(columns={"pfamA_acc":"template_domain"}).fillna("ND").drop_duplicates("clus_rep"), how="left")
    clus_df["size"] = clus_df["members"].apply(lambda x: len(x.split(",")))
    clus_df = clus_df.sort_values("size", ascending=False)
    clus_df.loc["same domain"] = False
    # Check if templates come from the same domain at least in 80% of the cases
    clus_df.loc[clus_df["template_domain"]!="ND", "same domain"] = clus_df.loc[clus_df["template_domain"]!="ND"]\
        .apply(lambda x: check_domain_consistency(x), 
        axis=1)

    # Save
    clus_df[["clus_rep", "template", "template_domain", "size", "same domain", "members"]].fillna("ND").to_csv("output/clusters.csv", index=False)


def check_domain_consistency(x):
    flag = False
    num_common_domains = len([i for i in str(x["members_template_domains"]).split(",") if i==x["template_domain"]])
    if num_common_domains > 0:
        fraction = num_common_domains / x["size"]
        if fraction > 0.8:
            flag = True
    return flag

def run_mmseqs_clustering(input_fasta, output_prefix, tmp_dir,
                          cluster_mode=2, cov_mode=0,
                          min_seq_id=0.25, coverage=0.9):
    command = [
        "mmseqs", "easy-cluster",
        "--cluster-mode", str(cluster_mode),
        "--cov-mode", str(cov_mode),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        input_fasta, output_prefix, tmp_dir
    ]
    
    subprocess.run(command, check=True)

if __name__ == "__main__":
    main()