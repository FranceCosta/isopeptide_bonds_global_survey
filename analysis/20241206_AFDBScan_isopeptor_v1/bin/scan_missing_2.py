#! /usr/env/python
# -*- coding: utf-8 -*-

"""

    Scan missing files (only those from split_0_214693 due to incorrect indexing, now fixed)
    Francesco Costa 2025-01-08 fcosta@ebi.ac.uk
    isopeptor version v0.0.75
    Run with:
    `conda activate isopeptor`
    `sbatch --job-name="isopeptor_afdb" -t 12:00:00 --mem 4GB -e log/%j.err -o log/%j.out --wrap="python3 bin/scan_missing_2.py"`

"""

import sys
sys.path.append("../../bin")
import os
import glob
from isopeptor.isopeptide import Isopeptide
import multiprocessing
import argparse
import time

DIR_LIST="/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v7/.files"

def main():
    
    with open(DIR_LIST, "rt") as fh:
        dirs = fh.readlines()
        output_directory = dirs[0].replace("\n", "")
    input_structures = os.path.join(output_directory, "cache_structures")
    jess(structures=input_structures, working_dir=output_directory)

def jess(structures:str, working_dir:str) -> None:
    """

        Runs Jess. Output in working_dir/jes_out.csv
        
        PARAMETERS
        ----------
        structures:str: list of structure file
        working_dir:str
    
        RETURNS
        -------
        Returns None
    
    """
    start = time.time()
    stdout_output = os.path.join(working_dir, f"jess_out.csv")
    i = Isopeptide("")
    # This is needed to pass structures as a list
    i.structure_files = [file.replace("\n", "") for file in open(structures, "rt").readlines()]
    i.predict()
    end = time.time()
    print(f"Writing output of {working_dir} after {int(end-start)} s. {len(i.structure_files)} files processed.")
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