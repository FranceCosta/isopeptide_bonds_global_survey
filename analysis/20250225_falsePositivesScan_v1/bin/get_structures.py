#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Download all structures associated with domains of
    interest
    Francesco Costa 2025-25-02 fcosta@ebi.ac.uk

"""


import sys
sys.path.append("../../bin")
from pfamenv import PFAM_HOST,PFAM_PASSWORD,PFAM_PORT,PFAM_USER,PFAM_VERSION
from mysql import connector
import numpy as np
import pandas as pd
import os
import requests

OUTPUT = "output"
DOMAINS = {"Helicase_C_2":"PF13307", 
            "Cu-oxidase_2":"PF07731"}

def main():

    # Get pfam info associated with each PDB
    columns = ["pdb_id", "pfamA_acc", "chain", "pdb_res_start", "pdb_res_end"]

    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()


    for pfamA_id in DOMAINS.keys():
        tmp_dir = os.path.join(OUTPUT, pfamA_id)
        os.makedirs(tmp_dir, exist_ok=True)
        # Get PDB ids
        cursor.execute(f"SELECT pdb_id,pfamA_acc,chain,pdb_res_start,pdb_res_end \
                        FROM {PFAM_VERSION}.pdb_pfamA_reg \
                        WHERE pfamA_acc = '{DOMAINS[pfamA_id]}'")

        output = cursor.fetchall()
        tmp_df = pd.DataFrame(output, columns = columns)
        
        # Download
        for index, row in tmp_df.iterrows():
            pdb_id = row["pdb_id"]
            chain = row["chain"]
            seq_start = row["pdb_res_start"]
            seq_end = row["pdb_res_end"]
            download_pdb(pdb_id, chain, seq_start, seq_end, save_dir=tmp_dir)
   
def download_pdb(pdb_id, chain, seq_start, seq_end, save_dir="."):
    """Download a PDB file from the RCSB PDB database and extract only the specified chain within a sequence range.
    
    Args:
        pdb_id (str): The PDB ID of the structure to download.
        chain (str): The chain of interest to extract.
        seq_start (int): The starting residue sequence number.
        seq_end (int): The ending residue sequence number.
        save_dir (str, optional): Directory where the PDB file should be saved. Defaults to current directory.
    
    Returns:
        str: Path to the downloaded PDB file with only the specified chain and sequence range, or None if the download failed.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, stream=True)
    
    if response.status_code == 200:
        os.makedirs(save_dir, exist_ok=True)
        file_path = os.path.join(save_dir, f"{pdb_id}_chain_{chain}_{seq_start}-{seq_end}.pdb")
        
        with open(file_path, "w") as f:
            for line in response.text.splitlines():
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if line[21] == chain:
                        try:
                            res_num = int(line[22:26].strip())
                            if seq_start <= res_num <= seq_end:
                                f.write(line + "\n")
                        except ValueError:
                            continue
                elif line.startswith("TER") or line.startswith("END"):
                    f.write(line + "\n")
        
        print(f"Downloaded {pdb_id}.pdb and saved chain {chain} from {seq_start} to {seq_end} to {file_path}")
        return file_path
    else:
        print(f"Failed to download {pdb_id}.pdb. HTTP status code: {response.status_code}")
        return None


if __name__ == "__main__":
    main()