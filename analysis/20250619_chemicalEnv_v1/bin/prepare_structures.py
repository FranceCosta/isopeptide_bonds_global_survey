"""
Prepare structures for foldMason alignment.
sbatch --job-name=prep_structs --output=log/%x_%j.out --error=log/%x_%j.err --time=12:00:00 --mem=8G --wrap="python bin/prepare_structures.py"
"""

import os
import pandas as pd
from Bio import SeqIO
import json
import biotite.database.afdb as afdb
import glob
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io import load_structure, save_structure
import sys
from dotenv import load_dotenv
load_dotenv("../../.env")

AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_ISOPEP_DOMAINS = os.getenv("AFDB_ISOPEP_DOMAINS")

isopep = ['Collagen_bind', 'GramPos_pilinBB', 'AgI_II_C2', 'Antigen_C',
       'Sgo0707_N2', 'DUF7926', 'DUF7929', 'DUF7925', 'DUF11', 'DUF7507',
       'DUF7619', 'DUF7933', 'GBS104-like_Ig', 'DUF7927', 'DUF7617',
       'SpaA', 'Cna_B', 'FctA', 'DUF5979', 'GramPos_pilinD1', 'DUF7601',
       'SpaA_4', 'SpaA_2', 'SpaA_3', 'GramPos_pilinD3', 'SdrD_B']


# Cache dir and file
STRUC_DIR = "/hps/nobackup/agb/research/francesco/20250624_isopep_bond_domains"
CACHED_SEQ_CSV = "tmp/sequences.csv"

def main():
    
    downloaded_structures_dir = os.path.join(STRUC_DIR, "downloaded_structures")
    cropped_structures_dir = os.path.join(STRUC_DIR, "cropped_structures")
    
    os.makedirs(downloaded_structures_dir, exist_ok=True)
    os.makedirs(cropped_structures_dir, exist_ok=True)

    str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
    afdb_isopep_domains = pd.read_csv(AFDB_ISOPEP_DOMAINS)

    str_df["phylum"] = str_df["taxonomy"].apply(lambda x: str(x).split(";")[1].replace(".", "") if len(str(x).split(";")) >= 2 else None)

    mapping_dict = str_df[["pfamA_acc", "pfamA_id"]].drop_duplicates().set_index("pfamA_id").to_dict()['pfamA_acc']

    df = str_df[(str_df["probability"]>0.65)&(~str_df["taxonomy"].isna())&(str_df["pfamA_id"].isin(isopep))]

    df = pd.concat([df, afdb_isopep_domains[afdb_isopep_domains["pfamA_id"].isin(isopep)]])
    df = df.sort_values("probability").drop_duplicates(["uniprot_acc", "seq_start", "seq_end"], keep='first')
    df["probability"] = df["probability"].fillna(0)

    if not os.path.isfile(CACHED_SEQ_CSV):
        records = (
            (record.id, str(record.seq))
            for record in SeqIO.parse(FULL_AFDB_SEQUENCES, format='fasta')
            if record.id in target_ids
        )
        seq_df = pd.DataFrame(records, columns=["uniprot_acc", "sequence"])
        seq_df.to_csv(CACHED_SEQ_CSV, index=False)
    else:
        seq_df = pd.read_csv(CACHED_SEQ_CSV)
    df = pd.merge(df, seq_df)

    # Adjust sequence borders
    df["seq_len"] = df["sequence"].apply(len)
    df["new_seq_start"] = df["seq_start"] - 30
    df["new_seq_end"] = df["seq_end"] + 30
    df.loc[df["new_seq_start"]<1, "new_seq_start"] = 1
    mask = (df["new_seq_end"]>df["seq_len"])
    df.loc[mask, "new_seq_end"] = df.loc[mask, "seq_len"]
    df["domain_sequence"] = df.apply(lambda x: x["sequence"][int(x["new_seq_start"]):int(x["new_seq_end"])], axis=1)  

    df["seq_start"] = df["seq_start"].astype(int)
    df["seq_end"] = df["seq_end"].astype(int)
    
    paths_dict = df[["struct_file", "uniprot_acc"]].dropna().drop_duplicates().set_index("uniprot_acc").to_dict()["struct_file"]
    for uniprot_acc in df["uniprot_acc"].unique():
        # If the path is not existing, download the structure from afdb
        if paths_dict.get(uniprot_acc, None):
            tmp_struc_path = paths_dict[uniprot_acc]
        else:
            tmp_struc_path = download_alphafold_model(uniprot_acc, downloaded_structures_dir)
        tmp_df = df[df["uniprot_acc"]==uniprot_acc]
        for _, row in tmp_df.iterrows():
            seq_start = row["seq_start"]
            seq_end = row["seq_end"]
            new_seq_start = row["new_seq_start"]
            new_seq_end = row["new_seq_end"]
            pfamA_acc = row["pfamA_acc"]
            out_path = os.path.join(cropped_structures_dir, pfamA_acc, f"{uniprot_acc}_{seq_start}_{seq_end}.pdb")
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            extract_region(tmp_struc_path, out_path, new_seq_start, new_seq_end)

def extract_region(in_path: str, out_path: str, seq_start: int, seq_end: int):
    structure = load_structure(in_path)
    region = structure[(structure.res_id >= seq_start) & (structure.res_id <= seq_end)]
    save_structure(out_path, region)

def download_alphafold_model(uniprot_acc: str, out_dir: str) -> str:
    out_path = os.path.join(out_dir, f"{uniprot_acc}.pdb")
    afdb.fetch(uniprot_acc, format="pdb", target_path=out_dir)

    return out_path

if __name__ == "__main__":
    main()