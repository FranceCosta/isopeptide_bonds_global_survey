#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Download PDB files containing isopeptide bonds
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
OUTPUT_DIR = "output"

def main():
    structures_dir = os.path.join(OUTPUT_DIR, "Positive_control")
    # Store here structures with no manual fixes
    unfixed_structures_dir = os.path.join(OUTPUT_DIR, "Positive_control_unfixed")

    df = pd.read_csv(TABLE)
    # Remove comments
    df = df[df["Chain"].isna() == False]  

    os.makedirs(structures_dir, exist_ok=True)
    os.makedirs(unfixed_structures_dir, exist_ok=True)

    # Download PDBs
    print("Downloading PDB files")
    pdbs = df["PDB code"].unique()
    plist = PDBList()
    for pdb in pdbs:
        plist.retrieve_pdb_file(pdb_code=pdb, file_format="pdb", overwrite=False, pdir=structures_dir)

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

    print("Ending")

if __name__ == "__main__":
    main()