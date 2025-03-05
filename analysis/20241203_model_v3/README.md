# Create model

`sbatch -t 24:00:00 --mem 32GB -J cluster --mail-type=END --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/cluster_matches.py"`