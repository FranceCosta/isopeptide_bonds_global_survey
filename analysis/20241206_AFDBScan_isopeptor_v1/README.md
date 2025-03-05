# Scan the entire AFDB with isopeptor

- first try using: `conda activate isopeptor` & `sbatch -t 500:00:00 --mem 128G --mail-type=END --mail-user=fcosta@ebi.ac.uk -J Isopeptor_AFDB --cpus-per-task=48 --wrap="python3 bin/scan_afdb.py"`; this code runs in one single node, using 48 jobs per time. With an anverage compute time of 3.12 days per job, this is not enough to complete the whole scan in a reasonable amount of time (less than 20 days);
- for this reason I have developed a second script that runs a job array with ~500 cuncurrent jobs: this should complete in less than a week. Remaining jobs will be run with a third script. I have also changed the set of template to use for this purpose by updating isopeptor. This leads to slightly different results and I have therefore re-run on the entire set using: `sbatch sh/scan_afdb_2.sh` which runs ovenight. The first set of structures was initially skipped. I have therefore run it with `sbatch --job-name="isopeptor_afdb" -t 12:00:00 --mem 4GB -e log/%j.err -o log/%j.out --wrap="python3 bin/scan_missing_2.py"`
- Annotate the final table using: `sbatch --job-name="annotate_afdb_scan" -t 48:00:00 --mem 16GB -e log/annotate_%j.err \
            -o log/annotate_%j.out --mail-type=END \
            --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/annotate.py"`


## Notes
The following taxonomy structures are there in the taxonomy annotation:
- For eukaryotes: superkingdom; kingdom; phylum; subphylum; clade (0-n, variable); class; infraclass; order; subfamily; genus
- for Bacteria: superkingdom; phylum; class; order; family; genus
- for archaea: superkingdom; phylum; clade (0-n, variable); class; order; family; genus