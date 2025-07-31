## Workflow

- `python bin/find_disulfide_bonds.py` to detect disulfide bonds in isopeptide bond-containing proteins (note that the "struct_file" column in the "AFDB_JESS_SCAN_TABLE" should contain the path to each AFDB structure in ".pdb" format);
- `python bin/format_table.py` to map disulfide bonds to isopeptide bond-containing Pfam domains;