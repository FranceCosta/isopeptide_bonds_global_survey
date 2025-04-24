#!/bin/bash -eu

# code from https://www.wwpdb.org/ftp/pdb-ftp-sites

rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/pdb/ output/PDB