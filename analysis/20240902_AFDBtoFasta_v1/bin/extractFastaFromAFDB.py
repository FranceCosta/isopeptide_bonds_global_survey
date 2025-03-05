#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Get fasta sequences from the AFDB structure files in a optimised way
    Version 3, use multiprocessing
    Francesco Costa 2024-10-07 fcosta@ebi.ac.uk

"""

import os
import glob
from multiprocessing import Pool

OUTPUT = "output/AFDB.fa"
AFDB_STRUCT_DIR = "/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v3/AFDB/structures/**/cleaned_structures_tmp/*.pdb" 
TMP_DIR = "/hps/nobackup/agb/research/francesco/tmp/AFDBFasta"
PROCESSES = 48
JOB_SIZE = 200000

aa_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
    'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'
}

def main():
    
    pool = Pool(PROCESSES)
    tmp_structures = []
    job_size = 0
    num_jobs = 0
    for pdb_file in glob.iglob(AFDB_STRUCT_DIR, recursive=True):
        if job_size < JOB_SIZE:
            tmp_structures.append(pdb_file)
            job_size += 1
        else:
            tmp_output_file = os.path.join(TMP_DIR, f"job_{num_jobs}.fa")
            pool.apply_async(convert_multiple, args=(tmp_structures, tmp_output_file))
            num_jobs += 1
            tmp_structures = []
            job_size = 0
    
    # Do for last bit outside job size
    if tmp_structures:
        tmp_output_file = os.path.join(TMP_DIR, f"job_{num_jobs}.fa")
        pool.apply_async(convert_multiple, args=(tmp_structures, tmp_output_file))
    
    # close() prevents any more tasks from being submitted
    pool.close()
    # join() waits for the worker processes to exit
    pool.join()
    
    # Merge all into output file
    with open(OUTPUT, 'w') as out:
        # Use iterator
        for tmp_file in glob.iglob(os.path.join(TMP_DIR, "*.fa"), recursive=True):
            with open(tmp_file, "r") as fh:
                    out.write(fh.read())

def convert_multiple(pdb_files:list, outfile:str)->None:
    with open(outfile, "w") as fh:
        for pdb_file in pdb_files:
            fh.write(convert_pdb(pdb_file))

def convert_pdb(pdb_file:str, chain="A")->str:
    seq = []
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[13:15].strip() == "CA" and line[21] == chain:
                resid = line[17:20].strip()
                seq.append(aa_dict.get(resid, ''))
    return f">{os.path.basename(pdb_file).replace('.pdb', '')}\n{''.join(seq)}\n"

if __name__ == "__main__":
    main()
