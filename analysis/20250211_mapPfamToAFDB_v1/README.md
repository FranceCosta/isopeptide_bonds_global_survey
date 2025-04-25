# Map to the whole AFDB the domain annotations of IB-containing domains

- `sbatch --job-name="map_pfam_to_afdb" -t 32:00:00 --mem 6GB -e map.err -o map.out --mail-type=END --ntasks 1  --cpus-per-task=48 --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/map.py"`

**Note that this code relies on a local installation of the [Pfam database](https://www.ebi.ac.uk/interpro/entry/pfam/#table)**