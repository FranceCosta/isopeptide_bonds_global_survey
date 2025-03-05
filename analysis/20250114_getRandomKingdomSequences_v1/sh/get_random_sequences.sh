# Select random proteins from pfamseq
PFAM_CONFIG="/nfs/research/agb/research/francesco/projects/20240212_isopeptideBonds_v1/20240529_findWithJess_v1/.pfamconfig"
PFAM_VERSION="pfam_37_0"
NUMBER_SEQS=10000

# Bacteria
mysql --defaults-extra-file=$PFAM_CONFIG $PFAM_VERSION --quick -N -e \
"
    SELECT pfamseq_acc,length,taxonomy 
    FROM pfamseq
    WHERE taxonomy LIKE 'Bacteria;%'
    AND length <= 1280
    ORDER BY RAND()
    LIMIT $NUMBER_SEQS;
" > output/bacteria_pfmaseq_random.tsv

# Bacteria
mysql --defaults-extra-file=$PFAM_CONFIG $PFAM_VERSION --quick -N -e \
"
    SELECT pfamseq_acc,length,taxonomy 
    FROM pfamseq
    WHERE taxonomy LIKE 'Archaea;%'
    AND length <= 1280
    ORDER BY RAND()
    LIMIT $NUMBER_SEQS;
" > output/archaea_pfmaseq_random.tsv