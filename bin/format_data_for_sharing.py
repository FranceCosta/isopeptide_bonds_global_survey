#! /usr/env/python
# -*- coding: utf-8 -*-

import os
from dotenv import load_dotenv
import pandas as pd
import shutil
import gzip
import zipfile
import tarfile
import datetime

load_dotenv(".env")

# Data to share
AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
BFVD_JESS_SCAN_TABLE = os.getenv("BFVD_JESS_SCAN_TABLE")
FA_PRED = os.getenv("FA_PRED")
AFDB_DOMAINS = os.getenv("AFDB_DOMAINS")
AF2_POS_SET_TEMPLATES = os.getenv("AF2_POS_SET_TEMPLATES")
PDB_BIOCHEM = os.getenv("PDB_BIOCHEM")
AF2_TEMPLATES_BIOCHEM = os.getenv("AF2_TEMPLATES_BIOCHEM")
BACTERIA_RANDOM_SEQUENCES = os.getenv("BACTERIA_RANDOM_SEQUENCES")
ARCHAEA_RANDOM_SEQUENCES = os.getenv("ARCHAEA_RANDOM_SEQUENCES")

SHARE_DIR = "to_share"

def main():
    
    # AFDB scan
    tmp_df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
    del tmp_df["chain"]
    del tmp_df["is_domain"]
    del tmp_df["struct_file"]
    tmp_df.to_csv(os.path.join(SHARE_DIR, "afdb_scan_isopeptor.csv"), index=False)

    # AFDB scan all domains
    tmp_df = pd.read_csv(AFDB_DOMAINS)
    tmp_df.to_csv(os.path.join(SHARE_DIR, "afdb_scan_isopeptor_pfam_annotations.csv"), index=False)

    # BFVD scan
    tmp_df = pd.read_csv(BFVD_JESS_SCAN_TABLE)
    del tmp_df["chain"]
    del tmp_df["is_domain"]
    del tmp_df["struct_file"]
    tmp_df.to_csv(os.path.join(SHARE_DIR, "bfvd_scan_isopeptor.csv"), index=False)

    # FA_pred
    tmp_df = pd.read_csv(FA_PRED)
    tmp_df.to_csv(os.path.join(SHARE_DIR, "fa_pred.csv"), index=False)

    # AF2 models (retain only ranked_0.pdb file!)
    output_zip = os.path.join(SHARE_DIR, "alphafold2_models.zip")
    output_gz = f"{output_zip}.gz"

    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(AF2_POS_SET_TEMPLATES):
            if "ranked_0.pdb" in files:  # Only proceed if the file exists in the directory
                file_path = os.path.join(root, "ranked_0.pdb")
                arcname = os.path.relpath(file_path, AF2_POS_SET_TEMPLATES)
                zipf.write(file_path, arcname)

    # Compress the zip file using gzip
    with open(output_zip, "rb") as f_in, gzip.open(output_gz, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Remove the original zip file
    os.remove(output_zip)

    # Biochem analysis
    tmp_df = pd.read_csv(PDB_BIOCHEM)
    del tmp_df["Notes"]
    del tmp_df["Interchain"]
    del tmp_df["Is bonded"]
    del tmp_df["Assigned bond"]
    del tmp_df["structure_path"]
    del tmp_df["match_residues"]
    del tmp_df["aro-isopep_planes_angle"]
    del tmp_df["water_resid_NZ"]
    del tmp_df["water_distance_NZ"]
    del tmp_df["water_NZ_ASA"]
    del tmp_df["water_OD_ASA"]
    del tmp_df["lys_x3"]
    del tmp_df["lys_x4"]

    tmp_df.to_csv(os.path.join(SHARE_DIR, "pdb_isopeptide_models_biochemical_analysis.csv"), index=False)

    tmp_df = pd.read_csv(AF2_TEMPLATES_BIOCHEM)
    del tmp_df["Notes"]
    del tmp_df["Interchain"]
    del tmp_df["Is bonded"]
    del tmp_df["Assigned bond"]
    del tmp_df["structure_path"]
    del tmp_df["match_residues"]
    del tmp_df["pdb_structure_path"]
    del tmp_df["lys_x3"]
    del tmp_df["lys_x4"]
    del tmp_df["af_r2"]

    tmp_df.to_csv(os.path.join(SHARE_DIR, "af2_isopeptide_models_biochemical_analysis.csv"), index=False)

    # Sequences
    shutil.copyfile(BACTERIA_RANDOM_SEQUENCES, os.path.join(SHARE_DIR, "afdb_bacteria_random_sequences.csv"))
    shutil.copyfile(ARCHAEA_RANDOM_SEQUENCES, os.path.join(SHARE_DIR, "afdb_archaea_random_sequences.csv"))
    
    # Compress all data
    tar = tarfile.open(os.path.join(SHARE_DIR, f"{datetime.date.today()}_data.tar.gz"), "w:gz")  
    for root, dirs, files in os.walk(SHARE_DIR):  
        for file in files:  
            if file != "README.md":  
                tar.add(os.path.join(root, file), arcname=os.path.relpath(os.path.join(root, file), SHARE_DIR))  
    tar.close()
    

if __name__ == "__main__":
    main()