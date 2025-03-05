#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Predict FAL probability
    `conda activate mlfa`
    `sbatch --job-name="predict_fa" -t 48:00:00 --mem 16GB -e predict_fa.err \
            -o predict_faout --mail-type=END \
            --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/predict_fa.py"`
    Francesco Costa 2025-02-11 fcosta@ebi.ac.uk

"""

import pandas as pd
import subprocess
import os
from dotenv import load_dotenv
load_dotenv("../../.env")
AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
AFDB_SCAN_SEQS = os.getenv("AFDB_SCAN_SEQS")

def main():
    # Create sequence file
    df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
    seq_df = pd.read_csv(AFDB_SCAN_SEQS)
    uniprot_accs = df[df["probability"]>.65]["uniprot_acc"].unique()
    with open("output/sequences.fasta", "wt") as fh:
        for index, row in seq_df[seq_df["uniprot_acc"].isin(uniprot_accs)].iterrows():
            fh.write(">{}\n{}\n".format(row["uniprot_acc"], row["sequence"]))

    # Run FA prediction
    run_fal_prediction(
        fasta_seqs="output/sequences.fasta",
        treks_dir="/hps/software/users/agb/research/francesco/software/FAL_prediction/Colab/scripts/",
        iupred_dir="/nfs/research/agb/research/vivian/software_vivian/software/", # This needs to be Vivian's version
        analysis_folder="output/analysis",
        results_folder="output/results",
        jobname="idp"
    )

def run_fal_prediction(fasta_seqs, treks_dir, iupred_dir, analysis_folder, results_folder, jobname):
    """
    Runs the FAL_prediction script with the given parameters.

    :param fasta_seqs: Path to the FASTA sequences file.
    :param treks_dir: Path to the Treks script directory.
    :param iupred_dir: Path to the IUPred software directory.
    :param analysis_folder: Path to the analysis output folder.
    :param results_folder: Path to the results output folder.
    :param jobname: Name of the job.
    """
    command = [
        "python3.7", 
        "/hps/software/users/agb/research/francesco/software/FAL_prediction/ML_predict.py", "predict",
        "--fasta_seqs", fasta_seqs,
        "--treks_dir", treks_dir,
        "--iupred_dir", iupred_dir,
        "--analysisfolder", analysis_folder,
        "--resultsfolder", results_folder,
        "--jobname", jobname
    ]
    
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("FAL Prediction completed successfully.")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error during FAL Prediction:")
        print(e.stderr)

if __name__ == "__main__":
    main()