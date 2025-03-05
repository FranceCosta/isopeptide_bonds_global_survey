#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Calculate biochemical properties on PDB isopeptide bond structures
    Francesco Costa 2024-10-25 fcosta@ebi.ac.uk

"""

import pandas as pd
import os
import numpy as np
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
from dotenv import load_dotenv

from dotenv import load_dotenv
load_dotenv("../../.env")
DATA_TABLE = os.getenv("TABLE")
POSITIVE_CONTROL = os.getenv("POSITIVE_CONTROL")
OUTPUT_TABLE = "output/pdb_biochem_properties.csv"

# ring atoms
AROMATICS = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"], 
             "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"], 
             "TRP": ["CD2", "CE2", "CE3", "CZ3", "CH2", "CZ2"]}

def main():
    df = pd.read_csv(DATA_TABLE)
    df = df[df["Chain"].isna() == False]
    #path_df = pd.read_csv(PATHS_TABLE, names=["structure_path"]).drop_duplicates()
    #path_df["PDB code"] = path_df["structure_path"].apply(lambda x: x.split("/")[-1].split("_")[0])
    #df = pd.merge(df, path_df, how="left")
    df["structure_path"] = df.apply(lambda x: os.path.join(POSITIVE_CONTROL, x["PDB code"].lower()+"_"+x["Chain"]+".pdb"), axis=1)
    df["match_residues"] = df.apply(lambda x: "_".join( 
                        [str(i) for i in sorted([x["Position 1\r\n(Bond 1)"], x["Position 2\r\n(catalytic)"], x["Position 3\r\n(Bond 2)"]])]
                        ), axis=1)
    # Restrict to intramolecular isopeptide bond
    cond1 = (df["Is bonded"] == True)
    cond2 = (df["Interchain"] == False)
    # Consider bad rotamers for comparison
    #cond3 = (df["Bad rotamer"] == False)
    df = df[cond1 & cond2]

    # Get Aromatic data
    aro_df = get_aro_params(df)
    df = pd.merge(df, aro_df, on=["structure_path", "Chain", "Position 1\r\n(Bond 1)"], how="left")
    
    # Get water to Nz and OE1/OD1 in bond 
    w_df = water_nz(df)
    df = pd.merge(df, w_df, how="outer", on=["PDB code", "Chain", "Position 1\r\n(Bond 1)"])
    w_df = water_od(df)
    df = pd.merge(df, w_df, how="outer", on=["PDB code", "Chain", "Position 1\r\n(Bond 1)"])
    
    # get water ASA
    df["water_NZ_ASA"], df["water_OD_ASA"] = zip(*df.apply(lambda x: get_water_ASA(x), axis=1))

    # get closest non catalytic oxygen (within 4A)
    df["o_atom_ref"], df["o_atom_distance"] = zip(*df.apply(lambda x: closest_o(x), axis=1))

    # Get ASA
    df["rASA"] = df.apply(getASA, axis=1)

    # Get bond length
    df["bond_length"] = df.apply(getBondLength, axis=1)
    
    # Get torsion angles
    df["pseudo_omega"], df["pseudo_psi"], df["pseudo_phi"], df["lys_x3"], df["lys_x4"] = zip(*df.apply(angles, axis=1))

    # Cis/trans annot
    df["cis"] = False
    df.loc[(df["pseudo_omega"]>=-60)&(df["pseudo_omega"]<=60), "cis"] = True

    df.to_csv(OUTPUT_TABLE, index=False)

def get_aro_params(df, dist_threshold=10) -> pd.DataFrame:
    """

        Get data asscoiated to aro rings: lenght, angle between isopep and ring planes;
        Sort aro rings ascending by the distance to the NZ then iterate among them and output 
        either the first aro with a cap detected or the closer aro. 

        PARAMETERS
        ----------
        dist_threshold=10, distance threshold to detect aromatic residues nearby NZ lys

    """
    
    dist_threshold=10
    # Cutoffs for aro cap
    # Adopt some approximate values since I am not considering hydrogens
    # max pi-H distance + C-H distance
    #dist_cutoff = 4.3+1.09
    #angle_cutoffs = [70, 110]
    dist_cutoff = 6
    angle_cutoffs = [30, 150]

    outlist = []
    for index, row in df.iterrows():
        struct_path = row["structure_path"]
        chain = row["Chain"]
        r1 = row["Position 1\r\n(Bond 1)"]
        r3 = row["Position 3\r\n(Bond 2)"]
        pdb_file = pdb.PDBFile.read(struct_path)
        # Exclude water and consider the right chain
        structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.hetero==False and atom.chain_id == chain])
        # Get isopep atoms
        isopep_atoms = struc.array([atom for atom in structure if (atom.res_id == r1 and atom.atom_name == "NZ") or \
                                                              (atom.res_id == r3 and atom.atom_name in ["CG", "OD1", "OD2"])])
        assert len(isopep_atoms) == 3, f"Found the following atoms in isopep_atoms: {isopep_atoms}"
        # Atoms within dist_threshold A from lys NZ
        lys_nz = struc.array([atom for atom in isopep_atoms if atom.res_id==r1 and atom.atom_name=="NZ"])
        atoms_within_threshold = surroundingAtom(lys_nz, structure, dist_threshold)
        # Get aro residues
        aro_residues_within_threshold = list(set(struc.array([atom for atom in atoms_within_threshold if atom.res_name in AROMATICS]).res_id))
        aro_rings = struc.array([atom for atom in structure if atom.res_id in aro_residues_within_threshold \
                                and atom.atom_name in AROMATICS.get(atom.res_name, [])])
        # This is necessary to exclude incomplete residues which may not have an aro ring
        aro_residues_within_threshold = list(set(aro_rings.res_id))

        # Order aro rings by distance with ring
        def get_distance(res_id):
            """Get distance with aro ring"""
            aro_atoms = struc.array([atom for atom in aro_rings if atom.res_id==res_id])
            return np.mean(aro_atoms.distance)

        aro_residues_within_threshold = sorted(aro_residues_within_threshold, key=lambda res_id: get_distance(res_id))
        
        # Iterate until a cap is detected otherwise return the first aro residue
        aro_output = None
        for res_id in aro_residues_within_threshold:
            aro_cap = False
            aro_atoms = struc.array([atom for atom in aro_rings if atom.res_id==res_id])
            ring_center = struc.centroid(aro_atoms)
            nz = lys_nz.coord[0]
            # Calculate the distance between the N atom and the ring center
            distance = struc.distance(ring_center, nz)

            # Get angles between aro ring and isopep bond planes
            def normal_vector(points):
                # Create two vectors on the plane
                vec1 = np.subtract(points[1], points[0])
                vec2 = np.subtract(points[2], points[0])
                # Compute the cross product to find the normal vector
                return np.cross(vec1, vec2)

            # Compute normal vectors for both planes
            normal1 = normal_vector(isopep_atoms.coord)
            normal2 = normal_vector(aro_atoms.coord)

            # Normalize the normal vectors
            normal1 = normal1 / np.linalg.norm(normal1)
            normal2 = normal2 / np.linalg.norm(normal2)

            # Compute the dot product of the normal vectors
            dot_product = np.dot(normal1, normal2)
            
            # Clip the dot product to avoid numerical issues (values slightly beyond [-1, 1])
            dot_product = np.clip(dot_product, -1.0, 1.0)
            
            # Calculate the angle in radians and convert to degrees
            angle_radians = np.arccos(dot_product)
            angle_degrees = np.degrees(angle_radians)

            if distance < dist_cutoff and ((angle_degrees < angle_cutoffs[0]) or (angle_degrees > angle_cutoffs[1])):
                aro_cap = True
            # Set the first by default (the closest)
            if aro_output == None:
                aro_output = [struct_path, chain, r1, res_id, aro_atoms.res_name[0], distance, angle_degrees, aro_cap]
            # If aro cap is detected, overwrite and break
            if aro_cap:
                aro_output = [struct_path, chain, r1, res_id, aro_atoms.res_name[0], distance, angle_degrees, aro_cap]
                break
        outlist.append(aro_output)

    return pd.DataFrame(outlist, 
            columns=["structure_path", "Chain", "Position 1\r\n(Bond 1)", 
            "aro_res_id", "aro_res_name", "distance_to_aro", "aro-isopep_planes_angle", "aro_cap"])     

def isopep_aro_plane_angle(row) -> float:
    """

        Get angle between the aromatic ring plane and the isopeptide bond plane (ASN OD1, CG, LYS NZ)

    """
    struct_path = row["structure_path"]
    chain = row["Chain"]
    r1 = row["Position 1\r\n(Bond 1)"]
    r3 = row["Position 3\r\n(Bond 2)"]
    aro_res_id = row["aro_res_id"]
    aro_res_name = row["aro_res_name"]
    pdb_file = pdb.PDBFile.read(struct_path)
    structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.hetero==False and atom.chain_id == chain])

    isopep_atoms = struc.array([atom for atom in structure if (atom.res_id == r1 and atom.atom_name == "NZ") or \
                                                              (atom.res_id == r3 and atom.atom_name in ["CG", "OD1", "OD2"])])
    assert len(isopep_atoms) == 3, f"Found the following atoms in isopep_atoms: {isopep_atoms}"
    aro_atoms = struc.array([atom for atom in structure if atom.res_id == aro_res_id and \
                                                        atom.atom_name in AROMATICS[atom.res_name][:3]])

    assert len(isopep_atoms) == 3, f"Found the following atoms in aro_atoms: {isopep_atoms}"

    def normal_vector(points):
        # Create two vectors on the plane
        vec1 = np.subtract(points[1], points[0])
        vec2 = np.subtract(points[2], points[0])
        # Compute the cross product to find the normal vector
        return np.cross(vec1, vec2)

    # Compute normal vectors for both planes
    normal1 = normal_vector(isopep_atoms.coord)
    normal2 = normal_vector(aro_atoms.coord)

    # Normalize the normal vectors
    normal1 = normal1 / np.linalg.norm(normal1)
    normal2 = normal2 / np.linalg.norm(normal2)

    # Compute the dot product of the normal vectors
    dot_product = np.dot(normal1, normal2)
    
    # Clip the dot product to avoid numerical issues (values slightly beyond [-1, 1])
    dot_product = np.clip(dot_product, -1.0, 1.0)
    
     # Calculate the angle in radians and convert to degrees
    angle_radians = np.arccos(dot_product)
    angle_degrees = np.degrees(angle_radians)

    return angle_degrees

def closest_o(row, threshold=5):
    """

        Find the closest oxygen atom to O atom

    """
    struct_path = row["structure_path"]
    r1 = row["Position 1\r\n(Bond 1)"]
    r2 = row["Position 2\r\n(catalytic)"]
    r3 = row["Position 3\r\n(Bond 2)"]
    res_3_aa = row["Residue 3"]
    chain = row["Chain"]
    atom_array = pdb.PDBFile.read(struct_path).get_structure()[0]
    structure = struc.array([atom for atom in atom_array if atom.chain_id == chain and atom.hetero == False])
    source = structure[(structure.res_id == r3) & (structure.chain_id == chain) & ((structure.atom_name=="OD1") | (structure.atom_name == "OD2"))]
    o_atom_ref, o_atom_distance = [None]*2
    if len(source) == 1:
        target = structure[(structure.element == "O") & (structure.res_id != r2) & (structure.res_id != r3)]
        if len(target) > 0:
            res_below_treshold = surroundingAtom(source, target, threshold)
            if len(res_below_treshold) > 0:
                sublist = []
                for o in res_below_treshold:
                    sublist.append([f"{o.res_name}_{o.res_id}_{o.atom_name}", struc.distance(source, o)[0]])
                # Report only the closest oxygen
                sublist = sorted(sublist, key=lambda x: x[-1])
                o_atom_ref, o_atom_distance = sublist[0]

    return o_atom_ref, o_atom_distance

def water_nz(df, threshold = 5):
    """

        Get the closer water molecule within threshold from isopep bond nz
    
    """
    outlist = []
    for _, row in df[["PDB code", "Chain", "structure_path"]].drop_duplicates().iterrows():
        pdb_id = row["PDB code"]
        chain = row["Chain"]
        protein = f"{pdb_id}_{chain}"
        struct_path = row["structure_path"]
        try:
            pdb_file = pdb.PDBFile.read(struct_path)
        except FileNotFoundError:
            continue
        # Consider the right chain
        structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.chain_id == chain])
        # Iterate over each isopep bond
        for _, row2 in df[(df["PDB code"] == pdb_id)&(df["Chain"] == chain)].iterrows():
            r1 = row2["Position 1\r\n(Bond 1)"]
            r2 = row2["Position 2\r\n(catalytic)"]
            r3 = row2["Position 3\r\n(Bond 2)"]
            try:
                if not structure[structure.res_name=="HOH"]:
                    continue
                source = structure[(structure.res_id == r1)&(structure.atom_name == "NZ")]
                target = struc.array([atom for atom in structure if atom.res_name == "HOH"])
                res_below_treshold = surroundingAtom(source, target, threshold)
                sublist = []
                for hoh in res_below_treshold:
                    sublist.append([pdb_id, chain, r1, hoh.res_id, struc.distance(source, hoh)[0]])
                # Report only the closest water
                sublist = sorted(sublist, key=lambda x: x[-1])
                if sublist:
                    outlist.append(sublist[0])
                
            except ValueError as e:
                ""

    return pd.DataFrame(outlist, columns=["PDB code", "Chain", "Position 1\r\n(Bond 1)", "water_resid_NZ", "water_distance_NZ"])

def water_od(df, threshold = 5):
    """

        Get the closer water molecule within threshold from isopep bond od1/2
    
    """
    outlist = []
    for _, row in df[["PDB code", "Chain", "structure_path"]].drop_duplicates().iterrows():
        pdb_id = row["PDB code"]
        chain = row["Chain"]
        protein = f"{pdb_id}_{chain}"
        struct_path = row["structure_path"]
        try:
            pdb_file = pdb.PDBFile.read(struct_path)
        except FileNotFoundError:
            continue
        # Consider the right chain
        structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.chain_id == chain])
        # Iterate over each isopep bond
        for _, row2 in df[(df["PDB code"] == pdb_id)&(df["Chain"] == chain)].iterrows():
            r1 = row2["Position 1\r\n(Bond 1)"]
            r2 = row2["Position 2\r\n(catalytic)"]
            r3 = row2["Position 3\r\n(Bond 2)"]
            if not structure[structure.res_name=="HOH"]:
                continue
            source = structure[(structure.res_id == r3)&((structure.atom_name == "OD1")|(structure.atom_name == "OD2"))]
            #print(structure[(structure.res_id == r3)])
            target = struc.array([atom for atom in structure if atom.res_name == "HOH"])
            res_below_treshold = surroundingAtom(source, target, threshold)
            sublist = []
            for hoh in res_below_treshold:
                sublist.append([pdb_id, chain, r1, hoh.res_id, struc.distance(source, hoh)[0]])
            # Report only the closest water
            sublist = sorted(sublist, key=lambda x: x[-1])
            if sublist:
                outlist.append(sublist[0])

    return pd.DataFrame(outlist, columns=["PDB code", "Chain", "Position 1\r\n(Bond 1)", "water_resid_OD", "water_distance_OD"])

def angles(row):
    """

        Get torsion angles

    """

    struct_path = row["structure_path"]
    r1 = row["Position 1\r\n(Bond 1)"]
    r3 = row["Position 3\r\n(Bond 2)"]
    res_3_aa = row["Residue 3"]
    chain = row["Chain"]
    atom_array = pdb.PDBFile.read(struct_path).get_structure()[0]
    
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
    
    except:
        ""

    return [pseudo_omega*180/np.pi, pseudo_psi*180/np.pi, pseudo_phi*180/np.pi, lys_x3*180/np.pi, lys_x4*180/np.pi]
       
def getBondLength(row) -> float:
    """

        Get isopeptide bond length

    """
    struct_path = row["structure_path"]
    r1 = row["Position 1\r\n(Bond 1)"]
    r3 = row["Position 3\r\n(Bond 2)"]
    res_3_aa = row["Residue 3"]
    chain = row["Chain"]
    atom_array = pdb.PDBFile.read(struct_path).get_structure()[0]
    
    # Get angles of bond and dihedral angle
    try:

        lys_nz = [atom for atom in atom_array if atom.res_id == r1 and atom.chain_id == chain and atom.atom_name == "NZ"][0]
        if res_3_aa == "N" or res_3_aa == "D":
            # CG is the last
            c1 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CG"][0]
        if res_3_aa == "E":
            # CD is the last
            c1 = [atom for atom in atom_array if atom.res_id == r3 and atom.chain_id == chain and atom.atom_name == "CD"][0]

        # This corresponds to the peptide omega angle (torsion angle on bond CN)
        dist = struc.distance(lys_nz, c1)
    except IndexError:
        # This occurs when the Lys does not have an N
        ""
    return dist

def surroundingAtom(source_atom:np.array, structure_array:list, threshold:float) -> list:
    """

        Find structure atoms within trehold distance in amstrongs from given atom

        PARAMETERS
        ----------
        source_atom: np.array:
        structure_array:list
        treshold:float

        RETURN
        ------
        resi_below_treshold:list: list of residues below given treshold
        
    """
    distances = struc.distance(source_atom, structure_array)
    structure_array.add_annotation("distance", dtype=float)
    structure_array.distance = distances
    resi_below_threshold = [atom for atom in structure_array if atom.distance < threshold]
    if resi_below_threshold:
        resi_below_threshold = struc.array(resi_below_threshold)
    return resi_below_threshold

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
    r1, r2, r3 = [int(i) for i in row["match_residues"].split("_")]
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
        ""
    except KeyError:
        print(row)

    return rASA

def get_water_ASA(row) -> float:
    """
    
        Calculates rASA of isopep residues. Returns average value
    
    """

    r1_asa, r2_asa = [np.NaN]*2
    struct_path = row["structure_path"]
    r1 = row["water_resid_NZ"]
    r2 = row["water_resid_OD"]
    pdb_file = pdb.PDBFile.read(struct_path)
    # Exclude hydrogens
    structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.element != "H" and ((atom.hetero==False) or (atom.atom_name=="HOH"))])
    # Consider whole structure to calculate sASA
    structure_sasa = struc.sasa(structure, point_number=500)
    # Get indeces
    if not np.isnan(r1):
        r1_indx = [i for i, atom in enumerate(structure) if atom.res_id == r1]
        r1_asa = sum([structure_sasa[i] for i in r1_indx])
    if not np.isnan(r2):
        r2_indx = [i for i, atom in enumerate(structure) if atom.res_id == r2]   
        r2_asa = sum([structure_sasa[i] for i in r2_indx])

    return [r1_asa, r2_asa]

if __name__ == "__main__":
    main()