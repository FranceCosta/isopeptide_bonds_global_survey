#!/bin/bash
# Download negative control

input="input/20241204_RSCB_search.txt"
output="output/Negative_control_3"

# Download structures of interest
sh/batch_download.sh -f $input -o $output -p

# Gunzip and replace
gunzip $output/*.gz