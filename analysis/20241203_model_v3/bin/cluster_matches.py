#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Cluster isopeptide bonds/signatures
    The first member from each cluster is the representative 
    and should be the only one considered
    Francesco Costa 2024-12-03 fcosta@ebi.ac.uk

"""

from Bio import SeqIO
from Bio import BiopythonParserWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import sys
sys.path.append("../../bin")
from cluster import sequenceAlignBiotite, cluster
import pandas as pd
import os
from dotenv import load_dotenv
load_dotenv("../../.env")

TABLE = os.getenv("TABLE")

STRUCTURE_PATHS = "output/structures"
CACHE_TABLE = f"output/jess_scan.csv"
OUTPUT = "output/cluster.csv"

ID_THRESHOLD = .3

def main():
    # Upload Jess scan results
    df = pd.read_csv(TABLE)
    res_df = pd.read_csv(CACHE_TABLE)
    
    #
    df["id"] = df.apply(lambda x: f"{x['PDB code']}_{x['Chain']}_"+
                    "_".join( 
                        [str(i) for i in sorted(
                            [
                            x["Position 1\r\n(Bond 1)"],
                            x["Position 2\r\n(catalytic)"],
                            x["Position 3\r\n(Bond 2)"]
                            ]
                            )]
                            ), axis=1)

    cond1 = (df["Is bonded"] == True)
    cond2 = (df["Interchain"] == False)
    cond3 = (df["Bad rotamer"] == False)
    cond4 = (df["Unusual geometry"]==False)
    cond5 = (df["resolution"]<=2.5)
    df = df[cond1 & cond2 & cond3 & cond4 & cond5]
    df["true_positive"] = 1

    # Get structure paths
    c_df = res_df.sort_values("rmsd").drop_duplicates(["target", "target_residues"], keep="first")
    with open(STRUCTURE_PATHS, "r") as fh:
        s_df = pd.DataFrame([l.strip() for l in fh.readlines()], columns = ["structure_path"])
        s_df["target"] = s_df["structure_path"].apply(lambda x: x.split("/")[-1].replace(".pdb", ""))
    
    c_df = pd.merge(c_df.assign(target=lambda x: x["target"].str.lower()), 
         s_df.assign(target=lambda x: x["target"].str.lower()), 
                     how="left")
    
    # Reduce sequences considering PDB start
    c_df["isopep_sequence"] = c_df.apply(reduceSeq, axis=1)
    
    c_df["id"] = c_df.apply(lambda x: x["target"]+"_"+"_".join(
        [str(i) for i in sorted([int(_) for _ in x["target_residues"].split("_")])]
        ),
        axis=1
    )
    
    # Give priority to true positives and low RMSD true negatives
    c_df = pd.merge(c_df, df[["id", "true_positive"]], how="left")
    c_df = c_df.sort_values(["true_positive", "rmsd"], ascending=[False, True])

    # Consider all sequences
    targets = c_df[["id", "isopep_sequence"]].set_index("id").to_dict()["isopep_sequence"]
    print("Number of sequences:")
    print(len(targets))
    
    outlist = []
    # .copy() is needed because otherwise it is not possible to pop entries from a dictionary while iterating it
    for t1 in targets.copy().keys():
        print(t1)
        for t2 in targets.keys():
            if t1 == t2:
                continue
            seq1 = targets[t1]
            seq2 = targets[t2]
            # Only max identity is calculated
            identity = sequenceAlignBiotite(seq1, seq2)
            
            outlist.append([t1, t2, identity])
        
        # Remove already done ones
        targets.pop(t1)
    
    seq_df = pd.DataFrame(outlist, columns=["Query", "Target", "identity"])
    clus_df = cluster(seq_df, ID_THRESHOLD)
    clus_df.to_csv(OUTPUT, index=False)

def reduceSeq(row):
    """
    
        Crop sequence in +- 20 amino acids from isopep signature considering PDB start position 
    
    """
    chain = row["chain"]
    residues = [int(i) for i in row["target_residues"].split("_")]
    min_res = min(residues)
    max_res = max(residues)

    # Consider chain as well
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonParserWarning)
        warnings.simplefilter('ignore', PDBConstructionWarning)
        seqs = list(SeqIO.parse(row["structure_path"], "pdb-atom"))
        for seq in seqs:
            if seq.id.split(":")[-1] == chain:
                break
                
    pdb_start = seq.annotations["start"]
    seq_start = min_res - pdb_start - 20
    seq_end = max_res - pdb_start + 20
    if seq_start < 0:
        seq_start = 0

    sequence = str(seq.seq)[seq_start:seq_end]
    
    return sequence



if __name__ == "__main__":
    main()