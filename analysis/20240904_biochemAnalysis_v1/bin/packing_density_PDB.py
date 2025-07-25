#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Calculate packing density on PDB isopeptide bond structures
    and void volume with and without isopeptide bond residues
    Francesco Costa 2025-07-17 fcosta@ebi.ac.uk

"""
import subprocess
import pandas as pd
import os
import numpy as np
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
from dotenv import load_dotenv
load_dotenv("../../.env")
DATA_TABLE = os.getenv("TABLE")

PROTEIN_VOLUME = "/nfs/research/agb/research/francesco/software/ProteinVolume_1"
POSITIVE_CONTROL = os.getenv("POSITIVE_CONTROL")
TMP_DIR = "tmp/protein_volume"
DOMAIN_BOUNDARIES = {
    "CnaA-like": [-10, 5],
    "CnaB-like": [-5, 30]
}
def main():
    os.makedirs(TMP_DIR, exist_ok=True)
    df = pd.read_csv(DATA_TABLE)
    df = df[df["Chain"].isna() == False]
    df["structure_path"] = df.apply(lambda x: os.path.join(POSITIVE_CONTROL, x["PDB code"].lower()+"_"+x["Chain"]+".pdb"), axis=1)
    cond1 = (df["Is bonded"] == True)
    cond2 = (df["Interchain"] == False)
    cond3 = (df["Isopeptide type"]!="Mutant")
    df = df[cond1 & cond2 & cond3]

    for index, row in df.iterrows():
        pdb_code = row["PDB code"]
        r1 = row["Position 1\r\n(Bond 1)"]
        r2 = row["Position 2\r\n(catalytic)"]
        r3 = row["Position 3\r\n(Bond 2)"]
        chain  = row["Chain"]
        bond_type = row["Isopeptide type"]
        struc_path = row["structure_path"]
        # Create boundaries
        seq_start =  r1 + DOMAIN_BOUNDARIES[bond_type][0]
        seq_end =  r3 + DOMAIN_BOUNDARIES[bond_type][1]
        pdb_file = pdb.PDBFile.read(struc_path)
        structure = pdb_file.get_structure()[0]
        structure = structure[(structure.hetero==False) & (structure.chain_id == chain) & (structure.atom_name!="H")]
        structure = structure[(structure.res_id >= seq_start) & (structure.res_id < seq_end)]
        struc_path = os.path.join(TMP_DIR, f"{pdb_code}_{seq_start}_{seq_end}.pdb")
        file = pdb.PDBFile()
        file.set_structure(structure)
        file.write(struc_path)

        p_vol = protein_volume(struc_path)
        print(pdb_code)
        print(p_vol)

        # Remove isopep bond residues and repeat
        backbone_atoms = ['N', 'CA', 'C', 'O']
        mask = np.isin(structure.res_id, [r1, r2, r3]) & ~np.isin(structure.atom_name, backbone_atoms)
        structure = structure[~mask]
        struc_path = os.path.join(TMP_DIR, f"{pdb_code}_{seq_start}_{seq_end}_no_isopep.pdb")
        file = pdb.PDBFile()
        file.set_structure(structure)
        file.write(struc_path)

        p_vol = protein_volume(struc_path)
        print(p_vol)

def protein_volume(pdb_path:str) -> list:
    """
    
        Run ProteinVolume to calculate protein volume

        PARAMETERS
        ----------
        pdb_path: str: path to pdb structure

        RETURNS
        -------
        list: 'Total_Volume(A3)', 'Void_Volume', 'VDW_Volume',  'Packing_Density', 'Time(ms)'

    """
    cmd = f"java -XX:+UseCompressedOops -XX:+UseSerialGC -jar {PROTEIN_VOLUME}/ProteinVolume_1.3.jar {pdb_path} | tail -n -1"
    volume = subprocess.check_output(cmd, shell=True, text=True)

    return [eval(v) for v in volume.split()[1:]]

if __name__ == "__main__":
    main()