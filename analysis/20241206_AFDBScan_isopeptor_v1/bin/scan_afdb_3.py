#! /usr/env/python
# -*- coding: utf-8 -*-

"""

    Run Jess on the whole AFDB
    Francesco Costa 2025-01-07 fcosta@ebi.ac.uk
    isopeptor version v0.0.78
    Run with:
    `conda activate isopeptor`
    `./sh/scan_afdb_2.sh`

"""

import sys
sys.path.append("../../bin")
import os
import glob
from isopeptor.isopeptide import Isopeptide
import multiprocessing
import argparse
import time

# - 1 is needed for correct indexing. This was initially not present
# which lead to a directory not being scanned
TASK_ID = int(os.getenv("SLURM_ARRAY_TASK_ID")) - 1
DISTANCE = 3

parser = argparse.ArgumentParser()

parser.add_argument(
    "--dir_list",
    help="Files containing the list of subdirectories.", 
    type=str,
)

def main():
    
    args = parser.parse_args()
    print(f"Task id is: {TASK_ID}")
    with open(args.dir_list, "rt") as fh:
        dirs = fh.readlines()
        output_directory = dirs[TASK_ID].replace("\n", "")
    input_structures = os.path.join(output_directory, "cache_structures")
    jess(structures=input_structures, working_dir=output_directory)

def jess(structures:str, working_dir:str, distance=DISTANCE) -> None:
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
    i = Isopeptide("", distance=distance)
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










