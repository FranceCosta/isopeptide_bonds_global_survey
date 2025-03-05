#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Run AF2 predictions of isopeptide bond structures 
    to compare with PDB ones
    Run with and without templates
    Francesco Costa 2024-10-24 fcosta@ebi.ac.uk
    `conda activate alphafold_v2`

"""

import os
import pandas as pd
import warnings
from Bio import SeqIO
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

TABLE = "data/20240923 Isopeptide, ester and thioester bond examples - Isopeptides FC modified.csv"
PDB_DIR = "output/Positive_control"
#OUT_DIR = "/nfs/research/agb/research/francesco/projects/20240212_isopeptideBonds_v1/20240529_findWithJess_v1/output/AF2"
OUT_DIR = "/nfs/research/agb/research/francesco/projects/20240212_isopeptideBonds_v1/20240529_findWithJess_v1/output/AF2_templates"
AF = "/nfs/research/agb/research/francesco/software/alphafold_non_docker/alphafold-2.3.2/"
DB = "/nfs/research/agb/research/vivian/2021-07-29-AlphaFold/databases/"
LOGS_DIR = "logs"

def main():
    df = pd.read_csv(TABLE)
    for _, row in df[["PDB code", "Chain"]].drop_duplicates().iterrows():
        pdb = row["PDB code"]
        chain = row["Chain"]
        output_structure = os.path.join(PDB_DIR, f"{pdb}_{chain}.pdb")
        sequence = list(SeqIO.parse(output_structure, "pdb-atom"))[0]
        # Replace X amino acids with G
        sequence.seq = sequence.seq.replace("X","G")
        sequence.id =  f"{pdb}_{chain}"
        sequence.description = ""
        seq_path = os.path.join(OUT_DIR, "sequences", f"{pdb}_{chain}.fa")
        with open(seq_path, "w") as fh:
            SeqIO.write(sequence, fh, "fasta")
        cmd = f'sbatch -t 06:00:00 --mem 64GB --gres=gpu:a100:1 \
                -o {os.path.join(LOGS_DIR, f"af2_{pdb}_{chain}.log")} \
                -e {os.path.join(LOGS_DIR, f"af2_{pdb}_{chain}.err")} \
                -J "af2_{pdb}_{chain}" \
                --wrap="cd {AF}; ./run_alphafold.sh -d {DB} -f {seq_path} -o {OUT_DIR} -t 2024-10-10"' # 1940-01-01
        os.system(cmd)


if __name__ == "__main__":
    main()