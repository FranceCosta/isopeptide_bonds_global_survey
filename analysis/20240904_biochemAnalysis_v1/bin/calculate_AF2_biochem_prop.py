#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Calculate biochemical properties on AF2 isopeptide bond structures
    AF2 structures have been obtained by running AF2 prediction
    Consider bad rotamers as well
    Automatically run for case with and without templates
    Francesco Costa 2024-10-25 fcosta@ebi.ac.uk

"""

import pandas as pd
import os
import numpy as np
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import warnings
from Bio import SeqIO
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

from dotenv import load_dotenv
load_dotenv("../../.env")
DATA_TABLE = os.getenv("TABLE")
POSITIVE_CONTROL = os.getenv("POSITIVE_CONTROL")

# ring atoms
AROMATICS = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"], 
             "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"], 
             "TRP": ["CD2", "CE2", "CE3", "CZ3", "CH2", "CZ2", "NE1",]}

def main():
    for AF_DIR, OUTPUT_TABLE in ([
        [os.getenv("AF2_POS_SET"), "output/af2_biochem_properties.csv"],
        [os.getenv("AF2_POS_SET_TEMPLATES"), "output/af2_templates_biochem_properties.csv"]
        ]):
        df = pd.read_csv(DATA_TABLE)
        df = df[df["Chain"].isna() == False]
        # Restrict to intramolecular isopeptide bond (NO)
        cond1 = (df["Is bonded"] == True)
        cond2 = (df["Interchain"] == False)
        #cond3 = (df["Bad rotamer"] == False)
        df = df[cond1&cond2]

        df["match_residues"] = df.apply(lambda x: "_".join( 
                            [str(i) for i in sorted([x["Position 1\r\n(Bond 1)"], x["Position 2\r\n(catalytic)"], x["Position 3\r\n(Bond 2)"]])]
                            ), axis=1)
        # Get af2 path
        df["structure_path"] = df.apply(lambda x: os.path.join(AF_DIR, x["PDB code"]+"_"+x["Chain"], "ranked_0.pdb"), axis=1)
        # Get PDB path (needed for matching residue correspondence between AF2 models and PDB structures)
        df["pdb_structure_path"] = df.apply(lambda x: os.path.join(POSITIVE_CONTROL, x["PDB code"].lower()+"_"+x["Chain"]+".pdb"), axis=1)
        
        # Match residues
        df["r1_af"], df["r2_af"], df["r3_af"] = zip(*df.apply(mapResidues, axis=1))
        # Manually fix wrong one
        df.loc[(df["PDB code"]=="6m3y")&(df["r1_af"]==467), "af_r2"] = 533

        # Get ASA
        df["rASA"] = df.apply(getASA, axis=1)

        # Get bond length
        df["bond_length"] = df.apply(getBondLength, axis=1)
        
        # Get torsion angles
        df["pseudo_omega"], df["pseudo_psi"], df["pseudo_phi"], df["lys_x3"], df["lys_x4"] = zip(*df.apply(angles, axis=1))

        df.to_csv(OUTPUT_TABLE, index=False)

def angles(row):
    """

        Get torsion angles

    """

    struct_path = row["structure_path"]
    r1 = row["r1_af"]
    r3 = row["r3_af"]
    res_3_aa = row["Residue 3"]
    chain = "A"
    atom_array = pdb.PDBFile.read(struct_path).get_structure()[0]
    
    pseudo_omega, pseudo_psi, pseudo_phi, lys_x3, lys_x4 = [np.NaN]*5
    # Get angles of bond and dihedral angle
    try:
        #lys_cb = [atom for atom in atom_array if atom.res_id == r1 and atom.chain_id == chain and atom.atom_name == "CB"][0]
        lys_cg = [atom for atom in atom_array if atom.res_id == r1 and atom.chain_id == chain and atom.atom_name == "CG"][0]
        lys_cd = [atom for atom in atom_array if atom.res_id == r1 and atom.chain_id == chain and atom.atom_name == "CD"][0]
        lys_ce = [atom for atom in atom_array if atom.res_id == r1 and atom.chain_id == chain and atom.atom_name == "CE"][0]
        lys_nz = [atom for atom in atom_array if atom.res_id == r1 and atom.chain_id == chain and atom.atom_name == "NZ"][0]
        if res_3_aa == "N" or res_3_aa == "D":
            # CG is the last
            c1 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CG"][0]
            # CB
            c2 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CB"][0]
            # CA
            c3 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CA"][0]
        if res_3_aa == "E":
            # CD is the last
            c1 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CD"][0]
            # CG
            c2 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CG"][0]
            # CB
            c3 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CB"][0]

        # This corresponds to the peptide omega angle (torsion angle on bond CN)
        pseudo_omega = struc.dihedral(c2, c1, lys_nz, lys_ce)
        # This corresponds to psi (one carbon is used instead of a second N)
        pseudo_psi = struc.dihedral(c3, c2, c1, lys_nz)
        # This corresponds to phi (lys_cd should be bound to oxygen and then to Nitrogrn to form the next peptide bond)
        pseudo_phi = struc.dihedral(c1, lys_nz, lys_ce, lys_cd)
        # Calc lysine dihedrals (X3 and X4)
        lys_x3 = struc.dihedral(lys_nz, lys_ce, lys_cd, lys_cg)
        lys_x4 = struc.dihedral(c1, lys_nz, lys_ce, lys_cd)
        pseudo_omega, pseudo_psi, pseudo_phi, lys_x3, lys_x4 = \
            [pseudo_omega*180/np.pi, pseudo_psi*180/np.pi, pseudo_phi*180/np.pi, lys_x3*180/np.pi, lys_x4*180/np.pi]
    
    except Exception as e:
        ""

    return [pseudo_omega, pseudo_psi, pseudo_phi, lys_x3, lys_x4]

def mapResidues(row):
    """

        Map isopeptide residues between PDB numbering and AF numbering

    """
    pdb_struct_path = row["pdb_structure_path"]
    r1 = row["Position 1\r\n(Bond 1)"]
    r2 = row["Position 2\r\n(catalytic)"]
    r3 = row["Position 3\r\n(Bond 2)"]
    sequence = list(SeqIO.parse(pdb_struct_path, "pdb-atom"))[0]
    pdb_seq_start = sequence.annotations["start"]
    r1_af  = r1 - pdb_seq_start + 1
    r2_af  = r2 - pdb_seq_start + 1
    r3_af  = r3 - pdb_seq_start + 1
    return [r1_af, r2_af, r3_af]

def getBondLength(row) -> float:
    """

        Get isopeptide bond length

    """
    struct_path = row["structure_path"]
    r1 = row["r1_af"]
    r3 = row["r2_af"]
    res_3_aa = row["Residue 3"]
    atom_array = pdb.PDBFile.read(struct_path).get_structure()[0]
    dist = np.NaN
    # Get angles of bond and dihedral angle
    try:

        lys_nz = [atom for atom in atom_array if atom.res_id == r1 and atom.atom_name == "NZ"][0]
        if res_3_aa == "N" or res_3_aa == "D":
            # CG is the last
            c1 = [atom for atom in atom_array if atom.res_id == r3 and atom.atom_name == "CG"][0]
        if res_3_aa == "E":
            # CD is the last
            c1 = [atom for atom in atom_array if atom.res_id == r3 and atom.atom_name == "CD"][0]

        # This corresponds to the peptide omega angle (torsion angle on bond CN)
        dist = struc.distance(lys_nz, c1)
    except IndexError:
        # This occurs when the Lys does not have an N
        ""
    return dist

def getASA(row) -> float:
    """
    
        Calculates rASA of isopep residues. Returns average value
    
    """
    MAX_ASA = { "rost_sander": { "LYS": 205, "ASP": 163, "GLU": 194, "ASN": 157 ,
    "SER": 130, # add for mutant
    "ALA": 106 # add for mutant
    }}
    rASA = np.NaN
    struct_path = row["structure_path"]
    r1 = row["r1_af"]
    r2 = row["r2_af"]
    r3 = row["r3_af"]
    #struct_path = os.path.join(PDB_DIR, f"{protein}.pdb")
    pdb_file = pdb.PDBFile.read(struct_path)
    # Exclude water
    structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.hetero==False and atom.element != "H"])
    # Consider whole structure to calculate sASA
    try:
        structure_sasa = struc.sasa(structure, point_number=500)
    except KeyError:
        ""
    try:
        # Get indeces
        r1_indx = [i for i, atom in enumerate(structure) if atom.res_id == r1]
        r2_indx = [i for i, atom in enumerate(structure) if atom.res_id == r2]
        r3_indx = [i for i, atom in enumerate(structure) if atom.res_id == r3]
        
        r1_aa = structure[r1_indx[0]].res_name
        r2_aa = structure[r2_indx[0]].res_name
        r3_aa = structure[r3_indx[0]].res_name
        
        r1_asa = sum([structure_sasa[i] for i in r1_indx]) / MAX_ASA["rost_sander"][r1_aa]
        r2_asa = sum([structure_sasa[i] for i in r2_indx]) / MAX_ASA["rost_sander"][r2_aa]
        r3_asa = sum([structure_sasa[i] for i in r3_indx]) / MAX_ASA["rost_sander"][r3_aa]
                
        rASA = [r1_asa, r2_asa, r3_asa]
    except IndexError:
        print(row)
        ""
    except KeyError:
        print(row)

    return rASA

if __name__ == "__main__":
    main()