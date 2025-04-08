#! /bin/bash -eu

# Run job array to scan the AFDB files generated with scan_afdb.py
# Run with
# `conda activate isopeptor`
# `./sh/scan_afdb_2.sh`


MAX_CUNCURRENT=500
#MAX_CUNCURRENT=5
# Get files from here
INPUT_FILES_DIR=/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v6/results
OUTPUT_DIR=/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v8
NUM_FILES=$( ls $INPUT_FILES_DIR | wc -l )

# Copy files
if [ ! -d $OUTPUT_DIR ]; then
    mkdir -p $OUTPUT_DIR
    cp -r $INPUT_FILES_DIR $OUTPUT_DIR
    # Remove old files
    rm $OUTPUT_DIR/results/*/jess_out.csv
fi

ls -d -1 $OUTPUT_DIR/results/** > $OUTPUT_DIR/.files

# --file_index is the index on the file list which indicate the file that should be run
sbatch --job-name="isopeptor_afdb" --array=1-$NUM_FILES%$MAX_CUNCURRENT \
    -t 167:30:00 --mem 4GB \
    -e log/%j.err -o log/%j.out \
    --wrap="python3 bin/scan_afdb_3.py  --dir_list=$OUTPUT_DIR/.files"
