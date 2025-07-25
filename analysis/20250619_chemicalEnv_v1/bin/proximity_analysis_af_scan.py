"""
Collect data about amino acids proximal to isopeptide bonds
In this version use pyjess to scan isopeptide bond positions to make
sure they are "isopeptide bond-equivalent" in negative structures
Do not remove isopep bond residues as they are part of the environment!

conda activate isopeptor
sbatch --time=12:00:00 --mem=16G --output=log/proximity_analysis_v2.log --error=log/proximity_analysis_v2.err --wrap="python bin/proximity_analysis_v2.py"

"""

import os
import numpy as np
import pandas as pd
from Bio.SeqUtils import seq3
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
from collections import Counter
from dotenv import load_dotenv
import pyjess
import glob

load_dotenv("../../.env")

AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_SCAN_SEQS = os.getenv("AFDB_SCAN_SEQS")
CLUSTER_REPS = "output/cluster_30.csv"
OUTPUT_CSV = "output/proximal_aa_afdb_clus_30.csv"
PYJESS_OUTPUT_CSV = "output/pyjess_scan.csv"
NEG_CSV = "output/no_isopep_mapped_to_hmm.csv"
NEG_STRUCTURES = "/hps/nobackup/agb/research/francesco/20250624_isopep_bond_domains/downloaded_structures"
TEMPLATES_DIR = "input/N_C_O_CA_templates/"
THRESHOLD = 6

isopep = ['Collagen_bind', 'GramPos_pilinBB', 'AgI_II_C2', 'Antigen_C',
       'Sgo0707_N2', 'DUF7926', 'DUF7929', 'DUF7925', 'DUF11', 'DUF7507',
       'DUF7619', 'DUF7933', 'GBS104-like_Ig', 'DUF7927', 'DUF7617',
       'SpaA', 'Cna_B', 'FctA', 'DUF5979', 'GramPos_pilinD1', 'DUF7601',
       'SpaA_4', 'SpaA_2', 'SpaA_3', 'GramPos_pilinD3', 'SdrD_B']

def main():
    str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
    clus_df = pd.read_csv(CLUSTER_REPS)

    neg_df = pd.read_csv(NEG_CSV) \
        [["pfamA_acc","uniprot_acc","seq_start", "seq_end","seq_pos1","seq_pos2","seq_pos3"]].rename(columns={
        "seq_pos1": "r1_bond",
        "seq_pos2": "r_cat",
        "seq_pos3": "r2_bond"
    })

    df = str_df[(str_df["probability"]>0.65)&(~str_df["taxonomy"].isna())&(str_df["pfamA_id"].isin(isopep))]

    # Add negative entries
    df = pd.concat([df, neg_df])
    df = df.sort_values("probability").drop_duplicates(["uniprot_acc", "seq_start", "seq_end"], keep='first')
    df["probability"] = df["probability"].fillna(0)

    df["id"] = df.apply(lambda row: f"{row.uniprot_acc}_{int(row.seq_start)}_{int(row.seq_end)}", axis=1)

    # Add missing info
    # Some structure paths are already present
    struct_file_dict = df[["struct_file", "uniprot_acc"]].drop_duplicates().dropna().set_index("uniprot_acc").to_dict()["struct_file"]
    df["struct_file"] = df["uniprot_acc"].map(struct_file_dict)
    # Other are not
    mask = df["struct_file"].isna()
    df.loc[mask, "struct_file"] = df.loc[mask].apply(lambda x: os.path.join(NEG_STRUCTURES, f'{x["uniprot_acc"]}.pdb'), axis=1)
    df["chain"] = "A"
    
    # Keep only cluster representatives
    df = pd.merge(df, clus_df, how="inner")
    df["r1_bond"] = df["r1_bond"].astype(int)
    df["r_cat"] = df["r_cat"].astype(int)
    df["r2_bond"] = df["r2_bond"].astype(int)

    # Run proximity analysis
    data = []
    for index, row in df.iterrows():
        struct_path = row["struct_file"]
        uniprot_acc = row["uniprot_acc"]
        id_ = row["id"]
        r1 = row["r1_bond"]
        r2 = row["r_cat"]
        r3 = row["r2_bond"]
        pdb_file = pdb.PDBFile.read(struct_path)
        structure = pdb_file.get_structure()[0]
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
        # Remove isopeptide residues-> no! Let's consider them for now
        #resi_below_threshold = resi_below_threshold[~np.isin(resi_below_threshold.res_id, [r1, r2, r3])]
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
        
        data.append([id_, uniprot_acc, r1, r2, r3]+list(aa_counter.values()))

    ch_df = pd.DataFrame(data, columns=["id", "uniprot_acc",  "r1_bond", "r_cat", "r2_bond"]+list(aa_counter.keys()))

    # Now collect pyjess data
    templates = [pyjess.Template.load(t, id=os.path.basename(t).replace(".pdb", "")) for t in glob.glob(TEMPLATES_DIR+"/*.pdb")]
    jess = pyjess.Jess(templates)
    distance = 2
    outlist = []
    for index, row in df[df["probability"]==0].drop_duplicates("uniprot_acc").iterrows():
        uniprot_acc = row["uniprot_acc"]
        struct_file = row["struct_file"]
        mol = pyjess.Molecule.load(struct_file, id=uniprot_acc)
        query = jess.query(mol, rmsd_threshold=distance, distance_cutoff=distance, max_dynamic_distance=distance, 
                        max_candidates=1000)
        for hit in query:
            atoms = hit.atoms()
            # This keeps the order of residues from templates (i.e.: bonded1, cat, bonded2)
            residue_numbers, residue_names = [], []
            for atom in atoms:
                if atom.residue_number not in residue_numbers:
                    residue_numbers.append(atom.residue_number)
                    residue_names.append(atom.residue_name)

            outlist.append([uniprot_acc, round(hit.rmsd, 3), int(residue_numbers[0]), int(residue_numbers[1]), int(residue_numbers[2])])
        
    j_df = pd.DataFrame(outlist, columns=["uniprot_acc", "rmsd", "r1_bond", "r_cat", "r2_bond"]) \
        .sort_values("rmsd", ascending=True).drop_duplicates(["uniprot_acc", "r1_bond", "r_cat", "r2_bond"], keep="first")
    j_df.to_csv(PYJESS_OUTPUT_CSV, index=False)
    ch_df = pd.merge(ch_df, j_df, how="left", on=["uniprot_acc", "r1_bond", "r_cat", "r2_bond"])
    ch_df.to_csv(OUTPUT_CSV, index=False)

if __name__ == "__main__":
    main()