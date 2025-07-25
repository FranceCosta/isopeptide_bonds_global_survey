import pandas as pd
import os
import shutil
import subprocess
from dotenv import load_dotenv
from pathlib import Path

load_dotenv("../../.env")

AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_SCAN_SEQS = os.getenv("AFDB_SCAN_SEQS")
TMP_DIR = Path("tmp/cluster")
shutil.rmtree(TMP_DIR)
TMP_DIR.mkdir(exist_ok=True)
OUTPUT_CSV = Path("output/cluster_30.csv")
OUTPUT_CSV.parent.mkdir(exist_ok=True)

FASTA_PATH = TMP_DIR / "domain_seqs.fasta"
DB_PATH = TMP_DIR / "seqDB"
CLUST_PATH = TMP_DIR / "clusterDB"
REP_SEQ_PATH = TMP_DIR / "rep_seqDB"
TSV_PATH = TMP_DIR / "cluster.tsv"

isopep = ['Collagen_bind', 'GramPos_pilinBB', 'AgI_II_C2', 'Antigen_C',
       'Sgo0707_N2', 'DUF7926', 'DUF7929', 'DUF7925', 'DUF11', 'DUF7507',
       'DUF7619', 'DUF7933', 'GBS104-like_Ig', 'DUF7927', 'DUF7617',
       'SpaA', 'Cna_B', 'FctA', 'DUF5979', 'GramPos_pilinD1', 'DUF7601',
       'SpaA_4', 'SpaA_2', 'SpaA_3', 'GramPos_pilinD3', 'SdrD_B']

CACHED_SEQ_CSV = "tmp/sequences.csv"
NEG_CSV = "output/no_isopep_mapped_to_hmm.csv"
def main():

    str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
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

    df = pd.merge(df, seq_df)

    df["seq_len"] = df["sequence"].apply(len)
    """
    # Elongate boundaries by 10 amino acids where necessary
    df["adj_seq_start"] = df["seq_start"]
    df["adj_seq_end"] = df["seq_end"]
    
    mask = (df["seq_start"]>df["r1_bond"])
    df.loc[mask, "adj_seq_start"] = df.loc[mask, "r1_bond"] - 10

    mask = (df["seq_end"]<df["r2_bond"])
    df.loc[mask, "adj_seq_end"] = df.loc[mask, "r2_bond"] + 10
    
    mask = (df["bond_type"] == "CnaB-like")
    df.loc[mask, "adj_seq_start"] = df.loc[mask, "r1_bond"] - 10
    df.loc[mask, "adj_seq_end"] = df.loc[mask, "r2_bond"] + 5

    mask = (df["bond_type"] == "CnaA-like")
    df.loc[mask, "adj_seq_start"] = df.loc[mask, "r1_bond"] - 5
    df.loc[mask, "adj_seq_end"] = df.loc[mask, "r2_bond"] + 30

    df["adj_seq_start"] = df["adj_seq_start"].apply(lambda x: 1 if x<=0 else x)
    df["adj_seq_end"] = df.apply(lambda x: x["seq_len"] if x["adj_seq_end"]>x["seq_len"] else x["adj_seq_end"], axis=1)

    df["domain_seq"] = df.apply(lambda x: x["sequence"][int(x["adj_seq_start"]) - 1:int(x["adj_seq_end"])], axis=1)
    """
    # Let's exclude domains missing any of the positions
    df = df[(~df["r1_bond"].isna()) & (~df["r_cat"].isna()) & (~df["r2_bond"].isna())]
    
    # Add sequence defined as +/-20 aa from first and last positions
    df["domain_seq"] = df.apply(lambda x: get_sequence(x), axis=1)

    #print(df[["seq_start", "seq_end", "r1_bond", "r_cat", "r2_bond", "seq_len"]])

    write_fasta(df, FASTA_PATH)

    run_mmseqs_clustering(FASTA_PATH, DB_PATH, CLUST_PATH, TMP_DIR)

    reps_df = extract_representatives(TSV_PATH)
    reps_df.to_csv(OUTPUT_CSV, index=False)

def get_sequence(row)->str:
    """

        Get sequence between min-20 and max+20 isopep signature residues from structure
    
    """
    r1 = row["r1_bond"]
    r2 = row["r_cat"]
    r3 = row["r2_bond"]
    seq_len = row["seq_len"]
    sequence = row["sequence"]
    seq_start = min([r1, r2, r3])
    seq_end = max([r1, r2, r3])
    
    seq_start = seq_start - 20
    seq_end = seq_end + 20
    if seq_start <= 0:
        seq_start = 1
    if seq_end > seq_len:
        seq_end=seq_len

    return str(sequence)[int(seq_start)-1:int(seq_end)-1]

def write_fasta(df, path):
    with open(path, "w") as f:
        for i, row in df.iterrows():
            header = f"{row.uniprot_acc}_{int(row.seq_start)}_{int(row.seq_end)}"
            f.write(f">{header}\n{row.domain_seq}\n")

def run_mmseqs_clustering(input_fasta, db_path, cluster_path, tmp_dir):
    subprocess.run(["mmseqs", "createdb", str(input_fasta), str(db_path)], check=True)
    subprocess.run([
        "mmseqs", "cluster", str(db_path), str(cluster_path), str(tmp_dir),
        "--min-seq-id", "0.3",
        "-c", "0.8",
        "--cov-mode", "0"
    ], check=True)
    subprocess.run(["mmseqs", "createtsv", str(db_path), str(db_path), str(cluster_path), str(TSV_PATH)], check=True)

def extract_representatives(tsv_path):
    # Only cluster representatives (first column, unique)
    tsv = pd.read_csv(tsv_path, sep="\t", header=None, names=["rep", "member"])
    reps = tsv["rep"].drop_duplicates().tolist()
    rep_info = []
    for rep in reps:
        acc, start, end = rep.split("_")
        rep_info.append({
            "uniprot_acc": acc,
            "seq_start": int(float(start)),
            "seq_end": int(float(end))
        })
    return pd.DataFrame(rep_info)

if __name__ == "__main__":
    main()
