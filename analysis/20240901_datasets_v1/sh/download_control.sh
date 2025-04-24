#!/bin/bash
# Download negative conntrol

input="input/20240611_negative_control.csv"
cache_file="tmp/list.txt"
output="output/Negative_control"

# Create input file with comma-separated structures
cat $input | tail -n +2 | cut -d "," -f2 | sort -u | paste -s -d, - | tr -d '"' > $cache_file

# Download structures of interest
../../../20240223_scanPDB_v1/sh/batch_download.sh -f $cache_file -o $output -p

# Gunzip and replace
gunzip $output/*.gz