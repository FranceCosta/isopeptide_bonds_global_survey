# Use Isopeptor to scan the [BFVD](https://academic.oup.com/nar/article/53/D1/D340/7906834)

- `sbatch --job-name="isopeptor_bfvd" -t 24:00:00 --mem 8GB -e log/%j.err -o log/%j.out --wrap="python3 bin/scan.py"`