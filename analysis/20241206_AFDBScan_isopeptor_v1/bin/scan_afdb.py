#! /usr/env/python
# -*- coding: utf-8 -*-

"""

    Run Jess on the whole AFDB
    Using modified templates from Set3 and other correct templates from set3
    `conda activate isopeptor`
    Francesco Costa 2024-12-06 fcosta@ebi.ac.uk
    isopeptor version v0.0.69

"""
import sys
sys.path.append("../../bin")
import os
import glob
from isopeptor.isopeptide import Isopeptide
import multiprocessing

OUTPUT_DIR = "/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v6"
#OUTPUT_DIR = "./tmp"
N_WORKERS = 48
N_JOBS = 1000
#N_WORKERS = 24
STRUCT_DIR = "/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v3/AFDB/structures"

def main():

    # Init workspace
    os.makedirs(os.path.join(OUTPUT_DIR, "results"), exist_ok=True)
    
    # Run Jess
    print("Collecting structures")
    #structures = list(glob.glob("/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v3/AFDB/structures/cleaned_structures_tmp/*.pdb", recursive=True))
    structures = list(glob.glob(f"{STRUCT_DIR}/**/cleaned_structures_tmp/*.pdb", recursive=True))
    print("Running Jess")
    jess_parallel(structures=structures, n_workers=N_WORKERS, n_jobs=N_JOBS, output_dir=os.path.join(OUTPUT_DIR, "results"))

def jess_wrapper(args):
    tmp_substructures, tmp_dir = args
    return jess(tmp_substructures, tmp_dir)

def jess_parallel(structures:list, n_workers:int, n_jobs:int, output_dir:str) -> None:
    """

        PARAMETERS
        ----------
        structures:list: list of structures full file paths
        n_workers:int
        n_jobs:int, number of jobs to send to the queue of n_workers
        output_dir:str
        
    """
    dirs = []
    substructures = []
    # Create multiple structures files each to be input to one Jess instance
    step = round(len(structures) / n_jobs)
    for i in range(0, len(structures), step):
        name = f"split_{i}_{i+step}"
        tmp_dir = os.path.join(output_dir, name)
        os.makedirs(tmp_dir, exist_ok=True)
        tmp_substructure = structures[i:i+step]
        # Cache structures
        structures_file = os.path.join(tmp_dir, "cache_structures")
        with open(structures_file, "w") as fh:
            for structure in structures[i:i+step]:
                fh.write(f"{structure}\n")
        substructures.append(tmp_substructure)
        dirs.append(tmp_dir)
    
    pool = multiprocessing.Pool(processes=n_workers)
    args = [(tmp_substructures, tmp_dir) for tmp_substructures, tmp_dir in zip(substructures, dirs)]
    r = pool.map(jess_wrapper, args)

    return None


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
    
    stdout_output = os.path.join(working_dir, f"jess_out.csv")
    #Path(stdout_output).touch()
    i = Isopeptide("")
    # This is needed to pass structures as a list
    i.structure_files = structures 
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












