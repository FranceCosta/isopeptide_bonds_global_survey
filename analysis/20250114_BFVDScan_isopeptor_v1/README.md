# Use Isopeptor to scan the [BFVD](https://academic.oup.com/nar/article/53/D1/D340/7906834)

- Scan the BFVD with `sbatch --job-name="isopeptor_bfvd" -t 24:00:00 --mem 8GB -e log/%j.err -o log/%j.out --wrap="python3 bin/scan.py"`
- Annotate the resulting hits with `sbatch --job-name="annotate_afdb_scan" -t 00:30:00 --mem 8GB -e log/annotate_%j.err -o log/annotate_%j.out --wrap="python3 bin/annotate.py"` **Note that this script relies on a local installtion of the [Pfam database](https://www.ebi.ac.uk/interpro/entry/pfam/#table)**