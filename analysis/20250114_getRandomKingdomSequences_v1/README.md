## Extract 1000 random sequences from Archaea and Bacteria

`sbatch -t 12:00:00 --mem 2G --mail-type=END --mail-user=fcosta@ebi.ac.uk  -J randomAFDBSeqs --wrap="python bin/get_random_sequences.py"`
- Fetch ~10k sequences from bacteria and archaea with length <= 1280 aa;
- match to AFDB sequences, at least 1000 entries;