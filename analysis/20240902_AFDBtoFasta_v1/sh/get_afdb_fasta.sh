#!/bin/bash -eu

# Convert protein structures to fasta sequences


#OUTPUT="output/AFDB.fa"
OUTPUT="tmp/AFDB_sh.fa"
#AFDB_STRUCT_DIR="/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v3/AFDB/structures/**/cleaned_structures_tmp/"
AFDB_STRUCT_DIR="/nfs/research/agb/research/francesco/projects/20240212_isopeptideBonds_v1/20240529_findWithJess_v1/tmp/test_pdbs/CnaA-like/*.pdb"

# Convert PDB to fasta
convertpdb() {
    
    pdb=$1

    declare -A aa_dict=(
    [ALA]=A
    [ARG]=R
    [ASN]=N
    [ASP]=D
    [CYS]=C
    [GLU]=E
    [GLN]=Q
    [GLY]=G
    [HIS]=H
    [ILE]=I
    [LEU]=L
    [LYS]=K
    [MET]=M
    [PHE]=F
    [PRO]=P
    [SER]=S
    [THR]=T
    [TRP]=W
    [TYR]=Y
    [VAL]=V
    )

    chain="A"

    printf ">%s\n" "$( basename $pdb '.pdb' )"
    
    for resid in  $( cat $pdb | grep ATOM | grep CA | grep $chain | cut -c 18-20 ); do
        printf "%s" "${aa_dict[$resid]}"
    done
    
    printf "\n"

}

export -f convertpdb

#convertpdb "/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v3/AFDB/structures/A0A023IDB7.pdb"

# Empty the output
> $OUTPUT 

find $AFDB_STRUCT_DIR -mindepth 0 -maxdepth 1 -name "*.pdb" \
     -exec bash -c 'convertpdb "$1"' - {} >> $OUTPUT \;