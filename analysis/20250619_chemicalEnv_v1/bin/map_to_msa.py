"""
Consider +/- 30 aa from domain boundaries
"""

import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio import AlignIO
import shutil
import sys
from dotenv import load_dotenv
load_dotenv("../../.env")

load_dotenv("../.env")

AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_DOMAINS= os.getenv("AFDB_DOMAINS")
AFDB_ISOPEP_DOMAINS = os.getenv("AFDB_ISOPEP_DOMAINS")
FULL_AFDB_SEQUENCES = os.getenv("FULL_AFDB_SEQUENCES")

isopep = ['Collagen_bind', 'GramPos_pilinBB', 'AgI_II_C2', 'Antigen_C',
       'Sgo0707_N2', 'DUF7926', 'DUF7929', 'DUF7925', 'DUF11', 'DUF7507',
       'DUF7619', 'DUF7933', 'GBS104-like_Ig', 'DUF7927', 'DUF7617',
       'SpaA', 'Cna_B', 'FctA', 'DUF5979', 'GramPos_pilinD1', 'DUF7601',
       'SpaA_4', 'SpaA_2', 'SpaA_3', 'GramPos_pilinD3', 'SdrD_B']


OUTPUT_CSV = "output/no_isopep_mapped_to_hmm.csv"
MSA_DIR = "output/MSAs"
CACHED_SEQ_CSV = "tmp/sequences.csv"
def main():
    str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
    afdb_isopep_domains = pd.read_csv(AFDB_ISOPEP_DOMAINS)
    df = str_df[(str_df["probability"]>0.65)&(~str_df["taxonomy"].isna())&(str_df["pfamA_id"].isin(isopep))]
    df = pd.concat([df, afdb_isopep_domains[afdb_isopep_domains["pfamA_id"].isin(isopep)]])
    df = df.sort_values("probability").drop_duplicates(["uniprot_acc", "seq_start", "seq_end"], keep='first')
    df["probability"] = df["probability"].fillna(0)
    target_ids = set(df["uniprot_acc"].unique())
    
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

    #AFDB_SCAN_SEQS = os.getenv("AFDB_SCAN_SEQS") # This is temporary for debug!
    #seq_df = pd.read_csv(AFDB_SCAN_SEQS)

    df = pd.merge(df, seq_df)

    df["seq_start"] = df["seq_start"].astype(int)
    df["seq_end"] = df["seq_end"].astype(int)

    df["seq_len"] = df["sequence"].apply(len)
    df["new_seq_start"] = df["seq_start"] - 30
    df.loc[df["new_seq_start"]<1, "new_seq_start"] = 1

    # Now for each domain find the positions along the HMM that correspond to isopep bonds and map them to the sequence
    results = []
    for pfamA_acc in df[df["pfamA_id"].isin(isopep)]["pfamA_acc"].unique():
        print(f"Mapping {pfamA_acc}")
        align = AlignIO.read(os.path.join(MSA_DIR, f"{pfamA_acc}_aa.fa"), "fasta")

        # Get columns that contain isopeptide bond positions
        data = []
        for _, row in df[df["pfamA_acc"] == pfamA_acc].iterrows():
            uniprot_acc = row["uniprot_acc"]
            new_seq_start, seq_start, seq_end = row["new_seq_start"], row["seq_start"], row["seq_end"]
            r1, rc, r2 = row["r1_bond"], row["r_cat"], row["r2_bond"]
            seq_id = f"{uniprot_acc}_{seq_start}_{seq_end}"
            new_seq_start = new_seq_start - 1
            rel_positions = [r1 - new_seq_start, rc - new_seq_start, r2 - new_seq_start]
            
            try:
                rec = next(r for r in align if seq_id in r.id)
            except:
                continue

            aligned = rec.seq
            pos_map = {}
            count = 0
            for i, res in enumerate(aligned):
                if res != "-":
                    count += 1
                    if count in rel_positions:
                        pos_map[count] = i

            mapped = [[p, pos_map[p], aligned[pos_map[p]]] for p in rel_positions if p in pos_map]
            if len(mapped) == 3:
                data.append([uniprot_acc, seq_start, seq_end] + [mapped[i][j] for i in range(3) for j in range(2)])

        t_df = pd.DataFrame(data, columns=["uniprot_acc", "seq_start", "seq_end", "r1", "r1_col", "r2", "r2_col", "r3", "r3_col"])
        # This is because some domains are split by insertions or are not well defined in pfam
        if len(t_df) == 0:
            print(f"No entries mapped for {pfamA_acc}")
            continue
        r1_col = t_df["r1_col"].value_counts().idxmax()
        r2_col = t_df["r2_col"].value_counts().idxmax()
        r3_col = t_df["r3_col"].value_counts().idxmax()

        # Get amino acids that are found at those column positions in domains where the isopeptide bond is absent
        for _, row in df[(df["pfamA_acc"] == pfamA_acc) & (df["probability"] == 0)].iterrows():
            uniprot_acc = row["uniprot_acc"]
            new_seq_start, seq_start, seq_end = row["new_seq_start"], row["seq_start"], row["seq_end"]
            new_seq_start = new_seq_start - 1
            seq_id = f"{uniprot_acc}_{seq_start}_{seq_end}"
            try:
                rec = next(r for r in align if seq_id in r.id)
            except:
                continue

            aligned = rec.seq
            aa_positions = []
            seq_idx = 0
            for i, res in enumerate(aligned):
                if res != "-":
                    seq_idx += 1
                if i in [r1_col, r2_col, r3_col]:
                    aa_positions.append((i, res, seq_idx if res != "-" else None))
            if len(aa_positions) == 3:
                results.append([pfamA_acc, uniprot_acc, seq_start, seq_end] + \
                    [aa for _, aa, _ in aa_positions] + [pos+new_seq_start for _, _, pos in aa_positions if pos != None])
    
    results_df = pd.DataFrame(results, columns=["pfamA_acc", "uniprot_acc", "seq_start", "seq_end",
            "aa1", "aa2", "aa3", "seq_pos1", "seq_pos2", "seq_pos3"])

    results_df.to_csv(OUTPUT_CSV, index=False)

if __name__ == "__main__":
    main()
