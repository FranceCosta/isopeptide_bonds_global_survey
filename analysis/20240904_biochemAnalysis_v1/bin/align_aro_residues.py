"""
Align aromatic residues to plot distribution of NZ atoms around them
"""

import pandas as pd
import numpy as np
import os
from dotenv import load_dotenv
from biotite.structure.io import pdb
import biotite.structure as struc
import shutil
import sys
sys.path.append("../../bin")
from sequence import get_sequence

load_dotenv("../../.env")
PDB_BIOCHEM = os.getenv("PDB_BIOCHEM")
OUTPUT_DIR = "output/aligned_aro"

def main():

    shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    pdb_df = pd.read_csv(PDB_BIOCHEM)

    # Get isopep sequence (between bonded residues)
    pdb_df["isopep_sequence"] = pdb_df.apply(get_sequence, pdb=True, axis=1)

    cond1 = (pdb_df["Is bonded"])
    cond2 = (pdb_df["Isopeptide type"]!="Mutant")
    cond3 = (~pdb_df["Unusual geometry/chemistry"])
    cond4 = (pdb_df["Resolution"]<=2.5)
    cond5 = (~pdb_df["Bad rotamer"])
    aro_df = pdb_df[cond1 & cond2 & cond3 & cond4 & cond5]
        
    # Save aligned versions of each aro-bond

    templates = {
        "trp": None,
        "tyr": None, 
        "phe": None
    }
    structure_path, chain = None, None

    for index, row in aro_df.sort_values("structure_path").iterrows():
        # Upload new structure if necessary
        tmp_structure_path = row["structure_path"]
        tmp_chain = row["Chain"]
        if tmp_structure_path != structure_path or tmp_chain != chain:
            structure_path = tmp_structure_path
            chain = tmp_chain
            pdb_file = pdb.PDBFile.read(structure_path)
            structure = struc.array([atom for atom in pdb_file.get_structure()[0] if atom.hetero==False and atom.chain_id == chain and atom.element != "H"])
        
        # Get aromatic residue and bond resid
        aro_res_id = row["aro_res_id"]
        aro_res_name = row["aro_res_name"].lower()
        r1 = row["Position 1\r\n(Bond 1)"]
        aro_cap = row["aro_cap"]
        # Upload template if necessary
        if not templates[aro_res_name]:
            template = structure[structure.res_id == aro_res_id]
            # Center ring centroir on origin
            ring_centroid = np.mean(template[~np.isin(template.atom_name, ["C", "N", "O", "CA", "CB"])].coord, axis=0)
            template.coord -= ring_centroid
            templates[aro_res_name] = template
            # Save template
            template = templates[aro_res_name]
            tmp_output = os.path.join(OUTPUT_DIR, f"template_{aro_res_name}.pdb")
            out_file = pdb.PDBFile()
            out_file.set_structure(template)
            out_file.write(tmp_output)
        
        template = templates[aro_res_name]
        target = structure[(structure.res_id == aro_res_id) | \
                            ((structure.res_id == r1) & (structure.atom_name == "NZ"))]

        tmp_template = template[~np.isin(template.atom_name, ["C", "N", "O", "CA", "CB"])]
        tmp_target = target[~np.isin(target.atom_name, ["C", "N", "O", "CA", "CB", "NZ"])]
        
        # Superpose the aro to the template
        fitted_target_aro, transform = struc.superimpose(tmp_template, tmp_target)
        # Transform aro+lys nz based on the superposition
        target = transform.apply(target)

        # Check if wrong orientation:
        backbone_atoms = target[np.isin(target.atom_name, ["CA", "CB"])]
        template_backbone = template[np.isin(template.atom_name, ["CA", "CB"])]
        if len(backbone_atoms) == 2 and len(template_backbone) == 2:
            vec_target = backbone_atoms.coord[1] - backbone_atoms.coord[0]
            vec_template = template_backbone.coord[1] - template_backbone.coord[0]
            # If wrong orientation, rotate around aro ring
            if np.dot(vec_target, vec_template) < 0:

                aro_ring_atoms = target[np.isin(target.atom_name, ["CG", "CZ"])]
                rotation_vector = aro_ring_atoms.coord[1] - aro_ring_atoms.coord[0]
                # Rotate by 180 degrees
                target = struc.rotate_about_axis(target, rotation_vector, angle=np.pi)

        # Save before
        tmp_output = os.path.join(OUTPUT_DIR, f"{aro_res_name}_aro_{aro_res_id}_lys_{r1}.pdb")
        if aro_cap:
            tmp_output = os.path.join(OUTPUT_DIR, f"{aro_res_name}_cap_aro_{aro_res_id}_lys_{r1}.pdb")
        out_file = pdb.PDBFile()
        out_file.set_structure(target)
        out_file.write(tmp_output)

if __name__ == "__main__":
    main()