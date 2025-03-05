#! /usr/env/python
# -*- coding: utf-8 -*-

"""

    Run Jess on the whole BFVD
    Francesco Costa 2025-01-14 fcosta@ebi.ac.uk
    isopeptor version v0.0.75
    Run with:
    `conda activate isopeptor`
    `sbatch --job-name="isopeptor_bfvd" -t 24:00:00 --mem 8GB -e log/%j.err -o log/%j.out --wrap="python3 bin/scan.py"`

"""
import sys
sys.path.append("../../bin")
import os
import glob
from isopeptor.isopeptide import Isopeptide
import multiprocessing

OUTPUT_DIR = "/hps/nobackup/agb/research/francesco/tmp/jessBFVD_v1"
STRUCT_DIR = "/hps/nobackup/agb/interpro/mblum/bfvd/structures"

def main():

    # Init workspace
    os.makedirs(os.path.join(OUTPUT_DIR, "results"), exist_ok=True)
    
    # Run Jess
    print("Running Jess")
    jess(OUTPUT_DIR, STRUCT_DIR)


def jess(working_dir:str, struct_dir:str) -> None:
    """

        Runs Jess. Output in working_dir/jes_out.csv
        
        PARAMETERS
        ----------
        working_dir:str
        struct_dir:str
    
        RETURNS
        -------
        Returns None
    
    """
    
    stdout_output = os.path.join(working_dir, f"jess_out.csv")
    i = Isopeptide(struct_dir=struct_dir)
    # This is needed to pass structures as a list
    i.predict()
    print(f"Writing output of {working_dir}")
    bonds = i.isopeptide_bonds
    with open(stdout_output, "w") as fh:
            if len(bonds) > 0:
                for bond in bonds:
                    row = [
                        bond.struct_file, bond.protein_name, str(bond.probability), bond.chain, str(bond.r1_bond), 
                        str(bond.r_cat), str(bond.r2_bond), 
                        bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name, bond.bond_type,
                        str(bond.rmsd), str(bond.r_asa), bond.template
                    ]
                    fh.write(",".join(row)+"\n")
        
    return None

if __name__ == "__main__":
    main()












