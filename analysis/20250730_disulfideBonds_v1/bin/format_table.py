import pandas as pd
import os
import numpy as np
from datetime import datetime
import sys
sys.path.append("../../")
from bin.constants import ISOPEP_DOMAINS
from dotenv import load_dotenv

load_dotenv("../../.env")
AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
DISULFIDE_CSV = "output/20250731_disulfide_bonds.csv"
DATE = datetime.today().strftime('%Y%m%d')
DISULFIDE_CSV = "output/disulfide_bonds.csv"
OUTPUT_CSV = f"output/{DATE}_mapped_disulfide_bonds.csv"

def main():

    str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE, low_memory=False)
    str_df = str_df[str_df["probability"]>.65]
    dis_df = pd.read_csv(DISULFIDE_CSV)
    str_df["superkingdom"] = str_df["taxonomy"].fillna("").apply(lambda x: x.split(";")[0].replace(".", ""))
    str_df["phylum"] = str_df["taxonomy"].fillna("").apply(lambda x: x.split(";")[1].replace(".", "") if len(x.split(";")) >= 2 else np.NaN)

    sub_df = pd.merge(str_df[str_df["pfamA_id"].isin(ISOPEP_DOMAINS)], 
        dis_df[["uniprot_acc", "res1", "res2"]].rename(columns={"res1":"cys_1", "res2": "cys_2"}),
        on=["uniprot_acc"], how="left")

    cond1 = ((sub_df["cys_1"]>sub_df["seq_start"]) & (sub_df["cys_1"]<sub_df["seq_end"]))
    cond2 = ((sub_df["cys_2"]>sub_df["seq_start"]) & (sub_df["cys_2"]<sub_df["seq_end"]))
    sub_df = sub_df[cond1 & cond2]

    sub_df["cys_1"] = sub_df["cys_1"].astype(int)
    sub_df["cys_2"] = sub_df["cys_2"].astype(int)

    sub_df.to_csv(OUTPUT_CSV, index=False)

if __name__ == "__main__":
    main()