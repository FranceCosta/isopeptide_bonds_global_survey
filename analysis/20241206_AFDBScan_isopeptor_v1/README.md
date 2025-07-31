# Scan the entire AFDB with isopeptor

## Workflow
- Run Isopeptor across the AFDB database: `sbatch sh/scan_afdb.sh` (requires AFDB database to be downloaded);
- Annotate the final table using: `python3 bin/annotate.py` (requires Pfam database);