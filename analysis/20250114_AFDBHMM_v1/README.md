# To reproduce:

1. `./sh/hmm.sh` 
2. `sbatch --job-name="format_HMM" -t 48:00:00 --mem 16GB -e log/format_%j.err\
                 -o log/format_%j.out --mail-type=END\
                 --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/formatHMM.py"`
3. `sbatch --job-name="annotate_HMM" -t 48:00:00 --mem 16GB -e log/annotate_%j.err\
                 -o log/annotate_%j.out --mail-type=END\
                 --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/annotateHMM.py"`