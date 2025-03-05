#!/bin/bash
# Download negative control

input="input/20240612_rcsb_pdb_custom_report.csv"
cache_file="tmp/list.txt"
output="output/Negative_control_2"

# Create input file with comma-separated structures
cat $input | tail -n +2 | cut -d "," -f2 | sort -u | paste -s -d, - | tr -d '"' > $cache_file

# Download structures of interest
sh/batch_download.sh -f $cache_file -o $output -p

# Gunzip and replace
gunzip $output/*.gz