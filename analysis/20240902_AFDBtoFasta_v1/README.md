## Fasta sequnces extraction

### Benchmark
- `sbatch -t 200:00:00 --mem 2G --mail-type=END --mail-user=fcosta@ebi.ac.uk  -J extractAFDB --wrap="python bin/extractFastaFromAFDB.py"` # This extracted 20 mln sequences in 8 days (10%)
- `sbatch -t 12:00:00 --mem 2G --mail-type=END --mail-user=fcosta@ebi.ac.uk  -J extractAFDB --wrap="./sh/get_afdb_fasta.sh"`

- Version 3 of python script: `sbatch -t 200:00:00 --ntasks=1 --cpus-per-task=48 --mem 32G --mail-type=END --mail-user=fcosta@ebi.ac.uk  -J extractAFDB --wrap="python bin/extractFastaFromAFDB.py"`

On a set of 44 structures:
- `extractFastaFromAFDB.py` time: 0.017s
- `get_afdb_fasta.sh` time: 0.303s

If this is still too slow, parallelize the python script. Done V