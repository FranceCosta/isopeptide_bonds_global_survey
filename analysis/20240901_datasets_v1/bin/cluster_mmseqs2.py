import pandas as pd
import os
import numpy as np
import subprocess
from dotenv import load_dotenv
import sys
from pathlib import Path

sys.path.append("../../bin")
from sequence import get_sequence

load_dotenv("../../.env")
DATA_TABLE = os.getenv("TABLE")
POSITIVE_CONTROL = os.getenv("POSITIVE_CONTROL")

TMP_DIR = "tmp"
OUTPUT_CSV = "output/clusters_mmseqs.csv"
Path(TMP_DIR).mkdir(exist_ok=True, parents=True)
Path("output").mkdir(exist_ok=True, parents=True)

MMSEQS_BIN = "mmseqs"

def main():
    df = pd.read_csv(DATA_TABLE)
    df = df[df["Chain"].notna()]
    df["structure_path"] = df.apply(lambda x: os.path.join(
        POSITIVE_CONTROL, f"{x['PDB code'].lower()}_{x['Chain']}.pdb"
    ), axis=1)

    df["match_residues"] = df.apply(lambda x: "_".join(
        [str(i) for i in[
            x["Position 1\r\n(Bond 1)"], 
            x["Position 2\r\n(catalytic)"], 
            x["Position 3\r\n(Bond 2)"]
        ]]
    ), axis=1)

    df["id"] = df.apply(lambda x: x["PDB code"]+"_"+x["Chain"]+x["match_residues"], axis=1)

    df = df[(df["Is bonded"] == True) & 
            (df["Interchain"] == False) & 
            (df["Bad rotamer"] == False)]

    df["isopep_sequence"] = df.apply(lambda x: get_sequence(x), axis=1)

    # Write FASTA for MMseqs
    fasta_path = os.path.join(TMP_DIR, "sequences.fasta")
    with open(fasta_path, "w") as fasta_file:
        for i, row in df.iterrows():
            fasta_file.write(f">{row['id']}\n{row['isopep_sequence']}\n")

    # Run clustering
    thresholds = [80, 60, 40, 20]
    cluster_assignments = {}

    for threshold in thresholds:
        tsv_path = run_mmseqs_clustering(fasta_path, TMP_DIR, "cluster", threshold)

        # Map from cluster representative to members
        clu_df = pd.read_csv(tsv_path, sep="\t", header=None, names=["rep", "member"])
        member_to_rep = clu_df.set_index("member")["rep"].to_dict()

        cluster_col = f"clus_rep_{threshold}"
        df[cluster_col] = df["id"].map(lambda seq_id: member_to_rep.get(seq_id, seq_id))
        cluster_assignments[threshold] = df[cluster_col]

    # Final output: id + cluster memberships
    df_out = df[["id"] + [f"clus_rep_{thr}" for thr in thresholds]]
    df_out.to_csv(OUTPUT_CSV, index=False)
    print(f"Written clustering results to {OUTPUT_CSV}")

def run_mmseqs_clustering(sequences_fasta, tmp_dir, output_prefix, min_seq_id):
    db_path = os.path.join(tmp_dir, f"db_{min_seq_id}")
    res_path = os.path.join(tmp_dir, f"res_{min_seq_id}")
    tdir = os.path.join(tmp_dir, f"tmp_{min_seq_id}")
    Path(tdir).mkdir(exist_ok=True)

    subprocess.run([
        MMSEQS_BIN, "createdb", sequences_fasta, db_path
    ], check=True)

    subprocess.run([
        MMSEQS_BIN, "cluster", db_path, res_path, tdir,
        "--min-seq-id", str(min_seq_id / 100.0),
        "-c", "0.8"         # minimum sequence coverage
    ], check=True)

    tsv_path = os.path.join(tmp_dir, f"cluster_rep_{min_seq_id}.tsv")
    subprocess.run([
        MMSEQS_BIN, "createtsv", db_path, db_path, res_path, tsv_path
    ], check=True)

    return tsv_path

if __name__ == "__main__":
    main()
