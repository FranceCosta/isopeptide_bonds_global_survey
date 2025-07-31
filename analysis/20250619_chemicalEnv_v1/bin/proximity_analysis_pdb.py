import os
import numpy as np
import pandas as pd
from Bio.SeqUtils import seq3
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
from collections import Counter
from dotenv import load_dotenv
load_dotenv("../../.env")

DATA_TABLE = os.getenv("TABLE")
POSITIVE_CONTROL = os.getenv("POSITIVE_CONTROL")
THRESHOLD = 6
OUTPUT_CSV = "output/proximal_aa_pdb.csv"

def main():
    df = pd.read_csv(DATA_TABLE)
    df = df[df["Chain"].isna() == False]
    df["structure_path"] = df.apply(lambda x: os.path.join(POSITIVE_CONTROL, x["PDB code"].lower()+"_"+x["Chain"]+".pdb"), axis=1)
    cond1 = (df["Is bonded"] == True)
    cond2 = (df["Interchain"] == False)
    cond3 = (df["Bad rotamer"] == False)
    cond4 = (df["Unusual geometry/chemistry"] == False)
    cond5 = (df["Isopeptide type"]!="Mutant")
    df = df[cond1 & cond2 & cond3 & cond4 & cond5]

    data = []
    for index, row in df.iterrows():
        struct_path = row["structure_path"]
        pdb_code = row["PDB code"]
        chain = row["Chain"]
        r1 = row["Position 1\r\n(Bond 1)"]
        r2 = row["Position 2\r\n(catalytic)"]
        r3 = row["Position 3\r\n(Bond 2)"]     
        isopep_type = row["Isopeptide type"]   
        pdb_file = pdb.PDBFile.read(struct_path)
        structure = pdb_file.get_structure()[0]
        structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.hetero==False and atom.chain_id == chain])
        r1_bb = structure[(structure.res_id == r1) & (structure.atom_name == "CA")]
        r2_bb = structure[(structure.res_id == r2) & (structure.atom_name == "CA")]
        r3_bb = structure[(structure.res_id == r3) & (structure.atom_name == "CA")]
        isopep_bb = r1_bb+r2_bb+r3_bb
        centroid = struc.centroid(isopep_bb)
        # Create virtual atom
        virtual_atom = struc.Atom(
            chain_id="A",
            res_id=-1,
            ins_code="",
            res_name="VRT",
            atom_name="VRT",
            element="X",
            coord=centroid,
            occupancy=1.0,
            b_factor=0.0
        )

        # Get distances
        distances = struc.distance(virtual_atom, structure)
        structure.add_annotation("distance", dtype=float)
        structure.distance = distances
        # Filter by threshold
        resi_below_threshold = structure[structure.distance < THRESHOLD]
        # Remove isopeptide residues
        resi_below_threshold = resi_below_threshold[~np.isin(resi_below_threshold.res_id, [r1, r2, r3])]
        # Remove backbone atoms except for ALA (consider "CB" in its case)
        resi_below_threshold = resi_below_threshold[(~np.isin(resi_below_threshold.atom_name, ["N", "C", "O", "CB", "CA"])) & \
            (resi_below_threshold.res_name!="ALA")] + resi_below_threshold[(~np.isin(resi_below_threshold.atom_name, ["N", "C", "O", "CA"])) & \
            (resi_below_threshold.res_name=="ALA")]
        sorted_indices = np.argsort(resi_below_threshold.res_id)
        resi_below_threshold = resi_below_threshold[sorted_indices]
        seq = ""
        if len(resi_below_threshold) > 0:
            seq = struc.to_sequence(resi_below_threshold)[0][0]
        aa_counter = Counter({aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY"})
        for aa in seq:
            aa_counter[aa] += 1
        
        data.append([pdb_code, chain, r1, r2, r3, isopep_type]+list(aa_counter.values()))

    ch_df = pd.DataFrame(data, columns=["PDB code", "Chain", "Position 1\r\n(Bond 1)", "Position 2\r\n(catalytic)", 
        "Position 3\r\n(Bond 2)", "Isopeptide type"]+list(aa_counter.keys()))
    ch_df.to_csv(OUTPUT_CSV, index=False)

if __name__ == "__main__":
    main()