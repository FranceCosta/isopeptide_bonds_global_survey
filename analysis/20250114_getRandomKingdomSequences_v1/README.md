## Extract 1000 random sequences from Archaea and Bacteria

## Workflow
- Fetch ~10k sequences from bacteria and archaea with length <= 1280 aa and match them to the AFDB database `python bin/get_random_sequences.py` (requires Pfam database and file with AFDB in fasta format, code to produce it in `../20240902_AFDBtoFasta_v1`);
