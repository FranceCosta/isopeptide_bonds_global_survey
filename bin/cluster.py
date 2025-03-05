#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Cluster sequences based on similarity
    Francesco Costa 2024-05-29 fcosta@ebi.ac.uk

"""

from Bio import SeqIO
import pandas as pd
from biotite.sequence import ProteinSequence
from biotite.sequence.align import SubstitutionMatrix, align_optimal

SUB_MATRIX = SubstitutionMatrix.std_protein_matrix()

def sequenceAlignBiotite(target_seq: str, subject_seq: str) -> float:
    """

        Performs global sequence alignment with Needleman-Wunsch algorithm and returns identity percentage calculated on aligned sequences
        note that terminal_penalty=False is necessary to allow the opening of gaps at beginning or ending of sequences

        PARAMETERS
        ----------
        target: str, protein target sequence;
        subject: str, protein subject sequence;

        RETURNS
        -------
        identity: float, max identity (subj-target vs target-subj);

    """
    def cleanSequence(sequence: str) -> str:
        """

            Remove unconventional caracthers from sequence

        """
        return sequence.replace('X', '').replace('U', '').replace('B', '').replace('Z', '')
    # Remove unconventional amino acids
    target_seq = cleanSequence(target_seq)
    subject_seq = cleanSequence(subject_seq)

    target = ProteinSequence(target_seq)
    subject = ProteinSequence(subject_seq)
    # blosum62 is used as substitution matrix
    ali = align_optimal(target, subject, SUB_MATRIX,
                        gap_penalty=-5, terminal_penalty=False)
    firts_ali = ali[0]  # all returned allignment have the same score
    id = 0
    for aa1, aa2 in zip(*firts_ali.get_gapped_sequences()):
        if aa1 == aa2 and aa1 != '-':
            id += 1

    return max(id/len(subject_seq), id/len(target_seq))

def cluster(seq_df: pd.DataFrame, treshold: float) -> pd.DataFrame:
    """
    
        Picks protein from seq all vs all identity DataFrame and removes from the DataFrame all sequences 
        within treshold of identity from the selected protein. Iterates until input DataFrame is empty. Proceeds in an
        ordered way
        
        PARAMETERS
        ----------
        seq_df: biotite generated all vs all output dataframe. Must contain 'Query', 'Target', 'identity';
        treshold: identity treshold to select targets;
        list of list with selected targets and coverage: [[target1, 0]...[targetn, N]] where N is the number of clusters;
        The first protein from each cluster is the representative member of the cluster
        
    """
    
    output, clustered = [], []
    protein = seq_df['Query'].to_list()[0]
    reduced_df = seq_df.copy()
    clus_num = 0
    
    while len(reduced_df) > 0:
        
        # Extract cluster, consider bidirectional identity
        tmp_df = reduced_df[((reduced_df["Query"]==protein) | (reduced_df["Target"]==protein)) & (reduced_df["identity"] >= treshold)]
        # Cluster rep is the first for each cluster
        proteins_in_cluster = [protein] + list(tmp_df["Target"].unique()) + list(tmp_df["Query"].unique())
        for prot in proteins_in_cluster:
            clustered.append(prot)
            output.append([prot, clus_num])
        
        # Remove proteins already clustered
        # This leaves out eventual protein that remains combined with the last cluster proteins but below treshold
        reduced_df = reduced_df[(reduced_df["Query"].isin(proteins_in_cluster) == False) & \
                                (reduced_df["Target"].isin(proteins_in_cluster) == False)]
        if len(reduced_df) == 0:
            # This fixes it
            for prot in seq_df['Query'].unique():
                c = 1
                if prot not in clustered:
                    output.append([prot, clus_num+c])
                    c+=1
            continue        
        # Select new protein
        protein = reduced_df['Query'].to_list()[0]
        clus_num += 1
    
    # Drop duplicates (cluster rep is found duplicated)
    return pd.DataFrame(output, columns=["Query", "cluster"]).drop_duplicates(keep='first')
