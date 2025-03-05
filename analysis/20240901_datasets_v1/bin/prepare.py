#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Download PDB files containing isopeptide bonds, 
    create non redundant datset based on sequence identity;
    prioritize manually reassigned NZ-structures and structures with side chains fixed by Rob;
    SPlit negative set into training and test
    Francesco Costa 2024-05-29 fcosta@ebi.ac.uk

"""

import pandas as pd
from Bio.PDB import PDBList
import shutil
from Bio import SeqIO
import os
import sys 
sys.path.append("../../bin")
from dotenv import load_dotenv
from cluster import cluster
from cluster import sequenceAlignBiotite
from structural import PDBcropper
import random
import warnings 
from Bio import BiopythonParserWarning
from Bio import BiopythonDeprecationWarning
from Bio.PDB.PDBParser import PDBConstructionWarning
warnings.simplefilter("ignore", BiopythonParserWarning)
warnings.simplefilter("ignore", BiopythonDeprecationWarning)
warnings.simplefilter("ignore", PDBConstructionWarning)

load_dotenv("../../.env")
TABLE = os.getenv("TABLE")
FIXED_PDB_DIR = os.getenv("FIXED_PDB_DIR")
MANUAL_PDB_REASSIGNED = os.getenv("MANUAL_PDB_REASSIGNED")
CLUSTER_TRESHOLD = 0.3
OUTPUT_DIR = "output"
NEGATIVE_SET = "output/Negative_control_3"

def main():
    structures_dir = os.path.join(OUTPUT_DIR, "Positive_control")
    # Store here structures with no manual fixes
    unfixed_structures_dir = os.path.join(OUTPUT_DIR, "Positive_control_unfixed")

    df = pd.read_csv(TABLE)
    # Remove comments
    df = df[df["Chain"].isna() == False]  

    os.makedirs(structures_dir, exist_ok=True)
    os.makedirs(unfixed_structures_dir, exist_ok=True)

    # Split Negative set into training and test
    os.makedirs(NEGATIVE_SET+"_test", exist_ok=True)
    os.makedirs(NEGATIVE_SET+"_training", exist_ok=True)
    files = os.listdir(NEGATIVE_SET)
    random.Random(42).shuffle(files)
    test_n = round(len(files)*.25)
    for file in files[:test_n]:
        shutil.copyfile(os.path.join(NEGATIVE_SET, file), os.path.join(NEGATIVE_SET+"_test", file))

    for file in files[test_n:]:
        shutil.copyfile(os.path.join(NEGATIVE_SET, file), os.path.join(NEGATIVE_SET+"_training", file))

    # Download PDBs
    print("Downloading PDB files")
    pdbs = df["PDB code"].unique()
    plist = PDBList()
    for pdb in pdbs:
        plist.retrieve_pdb_file(pdb_code=pdb, file_format="pdb", overwrite=False, pdir=structures_dir)
    
    # Downaload AlphaFold db models
    #print("Downloading AF2 files")
    #for uniprot_acc in df["UniProt Id"].unique():
    #    downloadAF2(uniprot_acc, AF_DIR)

    # Subselect chains
    print("Clustering sequences")
    for _, row in df[["PDB code", "Chain"]].drop_duplicates().iterrows():
        pdb = row["PDB code"]
        chain = row["Chain"]
        input_structure = os.path.join(structures_dir, f"pdb{pdb}.ent")
        # Overwrite if manually reassigned
        if os.path.exists(os.path.join(MANUAL_PDB_REASSIGNED, f"{pdb}.pdb")):
            input_structure = os.path.join(MANUAL_PDB_REASSIGNED, f"{pdb}.pdb")
        output_structure = os.path.join(structures_dir, f"{pdb}_{chain}.pdb")
        PDBcropper(input_structure, output_structure, 1, 10000000, chain_id=chain)
    
    # Backup bad rotamer structures
    for _, row in df[(df["Fixed"] == True) | (df["Bad rotamer"] == True)][["PDB code", "Chain"]].drop_duplicates().iterrows():
        pdb = row["PDB code"]
        chain = row["Chain"]
        input_structure = os.path.join(structures_dir, f"{pdb}_{chain}.pdb")
        output_structure = os.path.join(unfixed_structures_dir, f"{pdb}_{chain}.pdb")
        shutil.copyfile(input_structure, output_structure)

    # overwrite manually fixed chains
    print("Overwriting fixed chains")
    for _, row in df[df["Fixed"] == True][["PDB code", "Chain"]].drop_duplicates().iterrows():
        pdb = row["PDB code"]
        chain = row["Chain"]
        input_structure = os.path.join(FIXED_PDB_DIR, f"{pdb}.pdb")
        output_structure = os.path.join(structures_dir, f"{pdb}_{chain}.pdb")
        PDBcropper(input_structure, output_structure, 1, 10000000, chain_id=chain)

    # Extract sequences from PDB
    # Consider sequence surrounding the isopep bonds (+- 20 aa)
    # Give priority to good templates
    print("Clustering")
    df = df.sort_values(["Is bonded", "Bad rotamer", "Interchain"], ascending=[False, True, True])
    sequences = []
    for _, row in df.iterrows():
        pdb = row["PDB code"]
        chain = row["Chain"]
        r1 = row["Position 1\r\n(Bond 1)"]
        r2 = row["Position 2\r\n(catalytic)"]
        r3 = row["Position 3\r\n(Bond 2)"]
        seq_start = min([r1, r2, r3])
        seq_end = max([r1, r2, r3])
        output_structure = os.path.join(structures_dir, f"{pdb}_{chain}.pdb")
        sequence = list(SeqIO.parse(output_structure, "pdb-atom"))[0]
        residues = "_".join(
            [
                str(i) for i in sorted([r1, r2, r3])
            ]
        )
        sequence.id =  f"{pdb}_{chain}_{residues}"
        sequence.description = ""
        # Adjust seq start and end based on pdb seq structure start and end
        pdb_start = sequence.annotations["start"]
        seq_start = seq_start - pdb_start - 20
        seq_end = seq_end - pdb_start + 20
        if seq_start < 0:
            seq_start = 0
        sequence.seq = str(sequence.seq)[seq_start:seq_end]
        sequences.append(sequence)
    with open(os.path.join(OUTPUT_DIR, "input_sequences.fa"), "w") as fh:
        SeqIO.write(sequences, format="fasta", handle=fh)

    # Cluster sequences
    outlist = []
    partners = []
    for s1 in sequences:
        for s2 in sequences:
            if s1.id == s2.id:
                continue
            if [s1.id,s2.id] in partners:
                continue
            # The maximum is obtained
            identity = sequenceAlignBiotite(str(s1.seq), str(s2.seq))
            
            outlist.append([s1.id, s2.id, identity])
            partners.append([s1.id, s2.id])

    seq_df = pd.DataFrame(outlist, columns=["Query", "Target", "identity"])

    tmp_output = os.path.join(OUTPUT_DIR, f"pos_set_seq_identity.csv")
    seq_df.to_csv(tmp_output, index=False)
    clus_df = cluster(seq_df, CLUSTER_TRESHOLD)
    # The first protein of each cluster is the representative element
    clus_df.to_csv(os.path.join(OUTPUT_DIR, "cluster.csv"), index=False)
    print("Ending")

if __name__ == "__main__":
    main()





























