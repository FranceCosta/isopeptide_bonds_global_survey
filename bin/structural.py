#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Protein structure helper functions for AFDB run
    Francesco Costa 2024-05-29 fcosta@ebi.ac.uk

"""

from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import Select
import requests
import os

def downloadAF2(uniprot_acc: str, outdir: str, version='v4') -> bool:
    """

        Download AlphaFold models from the AlphaFold database and parse them to
        keep the structure of the desired seqeunce.

        PARAMETERS
        ----------
        auniprot_acc:str,  Uniprot accession code;
        outdir: str, path to directory where to save PDB structure;
        version: str, alphafold database version (default: v4);

        RETURNS
        bool: False if no structure or a wrong structure has been obtained;

    """
    download_link = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_{}.pdb"
    # print(f'Downloading {uniprot_acc}')
    output = os.path.join(outdir, f'{uniprot_acc}.pdb')

    try:
        # print("before")
        r = requests.get(
            download_link.format(uniprot_acc, version))
        # print("after")
        # print(f"r.status: {r.status_code}")
        r.raise_for_status()

    except requests.exceptions.HTTPError as e:
        # print(f"Found exception: {e}")
        return 0

    # write it to output
    open(output, 'wb').write(r.content)

    return 1


def PDBToTemplate(path_to_structure:str, path_to_output_structure:str, residue_indexes:list, chain_id='A') -> None:
    """
        Saves PDB substructure given chain start_pos and end_pos

        PARAMETERS
        ----------
        path_to_structure: str, 
        path_to_output_structure: str, 
        residue_indexes:list 
        chain_id: str, (A)

        RETURNS
        -------

        None, saves template in path_to_output_structure

    """
    class ResSelect(Select):
        def accept_residue(self, res):
            if res.id[1] in residue_indexes and res.parent.id == chain_id:
                return 1
            else:
                return 0

    s = PDBParser().get_structure('temp', path_to_structure)
    io = PDBIO()
    io.set_structure(s)
    io.save(path_to_output_structure, ResSelect())
    return None

def PDBcropper(path_to_structure: str, path_to_output_structure: str, start_pos: int, end_pos: int, chain_id='A') -> None:
    """
        Saves PDB substructure given chain start_pos and end_pos

        PARAMETERS
        ----------
        path_to_structure: str, 
        path_to_output_structure: str, 
        start_pos: int,
        end_pos: int, 
        chain_id: str, (A)

        RETURNS
        -------

        None, saves structure in path_to_output_structure

    """
    class ResSelect(Select):
        def accept_residue(self, res):
            if res.id[1] >= start_pos and res.id[1] <= end_pos and res.parent.id == chain_id:
                return 1
            else:
                return 0

    s = PDBParser().get_structure('temp', path_to_structure)
    io = PDBIO()
    io.set_structure(s)
    io.save(path_to_output_structure, ResSelect())
    return None