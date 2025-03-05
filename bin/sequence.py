#! /usr/env/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import re
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

S = re.compile('LP.T[G|A|N|D]')
S2 = re.compile('NP.TG')
S3 = re.compile('LP.GA')
S4 = re.compile('LA.TG')
S5 = re.compile('NPQTN')
S6 = re.compile('IP.TG')

def find_sortase(seq:str):
    """
    
        Find sortase motifs within 50 amino acids from C-ter
    
    """
    sortase = False
    seq = seq[-50:]
    if S.findall(seq):
        sortase = True
    elif S2.findall(seq):
        sortase = True
    elif S3.findall(seq):
        sortase = True
    elif S4.findall(seq):
        sortase = True
    elif S5.findall(seq):
        sortase = True
    elif S6.findall(seq):
        sortase = True
    
    return sortase

def get_sequence(row, pdb = True)->str:
    """

        Get sequence between r1_bond and r2_bond from structure
    
    """
    structure_path = row["structure_path"]
    if pdb:
        r1 = row["Position 1\r\n(Bond 1)"]
        r3 = row["Position 3\r\n(Bond 2)"]
    else:
        r1 = row["r1_af"]
        r3 = row["r3_af"]
    seq_start = min([r1, r3])
    seq_end = max([r1, r3])
    sequence = list(SeqIO.parse(structure_path, "pdb-atom"))[0]
    # Adjust seq start and end based on pdb seq structure start and end
    pdb_start = sequence.annotations["start"]
    return str(sequence.seq)[seq_start:seq_end]

