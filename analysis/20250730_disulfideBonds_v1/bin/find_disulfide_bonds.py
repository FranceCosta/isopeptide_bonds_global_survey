#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from Bio.PDB import NeighborSearch
import pandas as pd
from Bio import PDB
from datetime import datetime
from dotenv import load_dotenv

load_dotenv("../../.env")
parser = PDB.PDBParser(QUIET=True)

AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_STRUCTURES = "/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v3/AFDB/structures/**/cleaned_structures_tmp/*.pdb"
TMP_DIR = "/hps/nobackup/agb/research/francesco/tmp/IBD_dis"
ATOMS_RADII = {'S': 1.8}
DATE = datetime.today().strftime('%Y%m%d')
OUTPUT_CSV = f"output/{DATE}_disulfide_bonds.csv"

def main():
    
    # 
    os.makedirs(TMP_DIR, exist_ok=True)

    # Isopeptor-based AFDB scan
    str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE, low_memory=False)
    str_df = str_df[str_df["probability"]>.65]
    proteins_to_consider = str_df["struct_file"].unique()

    data = []
    for structure_path in proteins_to_consider:
        data.extend(process_structure(structure_path))
    
    df = pd.DataFrame(data, 
        columns=["uniprot_acc", "resname1", "resname2", 
                "atom1_name", "atom2_name", "res1", "res2", 
                "radii_sum", "distance", "atom1_plddt", "atom2_plddt"])
    df.to_csv(OUTPUT_CSV, index=False)
    
def process_structure(structure_path: str) -> list:
    data = []

    structure_name = os.path.basename(structure_path).replace(".pdb", "")
    structure = parser.get_structure("test", structure_path)
    cys_sg_atoms = [
        atom for atom in structure.get_atoms()
        if atom.get_id() == "SG" and atom.get_parent().get_resname() == "CYS"
    ]
    if len(cys_sg_atoms) >= 2:
        ns = NeighborSearch(cys_sg_atoms)
        close_atoms = ns.search_all(2*max(ATOMS_RADII.values()))

        for atom1, atom2 in close_atoms:
            res1 = atom1.get_full_id()[3][1]
            res2 = atom2.get_full_id()[3][1]
            
            # Exclude intra-residue and neighboring residues
            if abs(res1 - res2) <= 1:
                continue

            atom1_name = atom1.name
            atom2_name = atom2.name

            radii_sum = ATOMS_RADII[atom1.get_id()[0]] + ATOMS_RADII[atom2.get_id()[0]]
            distance = atom1 - atom2

            # Exclude if distance is above radii sum - clash threshold
            if distance > (radii_sum):
                continue

            resname1, resname2 = atom1.get_parent().get_resname(), \
                                    atom2.get_parent().get_resname()

            data.append([structure_name, resname1, resname2, 
                                atom1_name, atom2_name, res1, res2, 
                                radii_sum, distance, atom1.bfactor, atom2.bfactor])
        
    return data

if __name__ == "__main__":
    main()