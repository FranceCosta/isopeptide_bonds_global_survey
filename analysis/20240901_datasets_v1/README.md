# Format intramolecular isopeptide bonds dataset and download PDB structures

## Workflow
- Manually download `PDB_bad_bonds_check_reassigned` and `manually_reassigned` Data from: Costa,F. et al. (2025) Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures. Bioinformatics Advances, vbaf049. Place them in `../../data`
- Download missing PDB structures with `python bin/prepare.py`;
- Annotate isopeptide bond-containing proteins with Pfam domains: `python bin/add_pfam_annot.py` (this requires a local Pfam database)