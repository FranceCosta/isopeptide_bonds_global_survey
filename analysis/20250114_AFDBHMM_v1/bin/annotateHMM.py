#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Annotate sequences with new domains created in PFAM using HMMSCAN
    Francesco Costa 2024-07-22 fcosta@ebi.ac.uk
    Sequences are cached 
    `sbatch --job-name="annotate_HMM" -t 48:00:00 --mem 16GB -e log/annotate_%j.err\
                 -o log/annotate_%j.out --mail-type=END\
                 --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/annotateHMM.py"`

"""

import time
import os
import subprocess
import sys
import pandas as pd
from pathlib import Path
sys.path.append("../../bin")
from pfamenv import PFAM_HOST, PFAM_PASSWORD, PFAM_PORT, PFAM_USER, PFAM_VERSION
from mysql import connector
import re
from datetime import datetime
import numpy as np
import glob

DATE = datetime.today().strftime('%Y%m%d')
HMMSCAN = "/hps/software/users/agb/research/francesco/software/hmmer-3.3.2/bin/hmmscan"
HMMPRESS = "/hps/software/users/agb/research/francesco/software/hmmer-3.3.2/bin/hmmpress"
INPUT_TABLE = "output/formatted_HMM.csv"
TMP_SEQS = "tmp/hmm_sequences.fa"
TMP_HMM_DIR = "tmp/HMMs"
HMMSCAN_OUTPUT = "tmp/hmscan.out"
# Specify here which domains to be used for the scan
DOMAINS = ["PF24346", "PF24514", "PF24517", "PF24547", "PF24558", "PF24593", "PF24595"]
OUTPUT_TABLE = f"output/{DATE}_HMMScan.csv"

def main():

    os.makedirs(TMP_HMM_DIR, exist_ok=True)
    # Remove old HMM file
    if os.path.isfile(os.path.join(TMP_HMM_DIR, "HMM")):
        os.remove(os.path.join(TMP_HMM_DIR, "HMM"))
    
    # Get sequences that hmmscan will use
    hmm_df = pd.read_csv(INPUT_TABLE, low_memory=False)
    if not os.path.exists(TMP_SEQS):
        seq_df = getSequences(hmm_df["uniprot_acc"].unique())
        with open(TMP_SEQS, "w") as fh:
            for _, row in seq_df.iterrows():
                fh.write(">"+row["uniprot_acc"]+"\n"+row["sequence"]+"\n")

    # Download HMMs and concat
    domain_data = {}
    for domain in DOMAINS:
        
        Path(os.path.join(TMP_HMM_DIR, "HMM")).touch()
        # Download if not already existing
        if not os.path.exists(os.path.join(TMP_HMM_DIR, domain)):
            cmd = f"cd {TMP_HMM_DIR}; pfco {domain}"
            r = subprocess.run(cmd, shell=True, text=True, capture_output=True)
            if r.returncode != 0:
                print("Error encountered")
                print(r.stderr)
                sys.exit()
        
        # Get gathering thresholds
        ga = re.compile(r"GA.+(\d+\d+.\d+\d+).+(\d+\d+.\d+\d+);")
        ac = re.compile(r"ID   (.+)")
        cl = re.compile(r"CL   (CL\d+)")
        clan_acc = np.NaN
        with open(os.path.join(TMP_HMM_DIR, domain, "DESC"), "rt") as desc:
            for line in desc:
                if line.startswith("ID"):
                    pfamA_id = ac.findall(line)[0]
                if line.startswith("GA"):
                    thresholds = ga.findall(line)
                if line.startswith("CL"):
                    #print(line, cl.findall(line))
                    clan_acc = cl.findall(line)[0]

        d_thr, s_thr = float(thresholds[0][0]), float(thresholds[0][1])
        domain_data[domain] = {"pfamA_id": pfamA_id, "domain_threshold": d_thr, 
                                        "sequence_threshold": s_thr, "clan_acc": clan_acc}
        
        # Concat to HMM file
        with open(os.path.join(TMP_HMM_DIR, "HMM"), "a") as outfile:
            outfile.write("\n")
            with open(os.path.join(TMP_HMM_DIR, domain, "HMM"), "rt") as hmm:
                for line in hmm:
                    if "NAME  SEED" in line:
                        line = f"NAME  {domain}\n"
                    outfile.write(line)
    # Init database
    for index_file in glob.glob(os.path.join(TMP_HMM_DIR, "HMM.h*")):
        os.remove(index_file)

    cmd = f"{HMMPRESS} {os.path.join(TMP_HMM_DIR, 'HMM')}"
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if r.returncode != 0:
        print("Error encountered during index creation")
        print(r.stderr)
        sys.exit()

    # Run HMMscan
    cmd = f"{HMMSCAN} --domtblout {HMMSCAN_OUTPUT} {os.path.join(TMP_HMM_DIR, 'HMM')} {TMP_SEQS}"
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if r.returncode != 0:
        print("Error encountered during scan")
        print(r.stderr)
        sys.exit()

    # Parse results
    # NOTE: This code avoids overwriting existing domains
    pf_df = pd.read_csv(HMMSCAN_OUTPUT, comment='#', 
                          names=["target name", "accession1", "tlen", "query name", "accession2", "qlen", "E-value", "full_seq_score", "full_seq_bias",
                                 "#", "of", "c-Evalue", "i-Evalue", "domain_score", "bias", "hmm_from",
                                 "hmm_to", "ali_from", "ali_to", "from", "to", "acc", "description of target"], sep=r'\s+')
    
    # consider a domain specific GA
    p_df = pd.DataFrame()
    for pfamA_acc in domain_data:
        p_df = pd.concat([p_df, 
                          pf_df[(pf_df["full_seq_score"]>domain_data[pfamA_acc]["sequence_threshold"]) &\
                     (pf_df["domain_score"]>domain_data[pfamA_acc]["domain_threshold"]) &\
                     (pf_df["target name"] == pfamA_acc)]\
              [["target name", "query name", "from", "to"]]\
                   .rename(columns={"target name": "pfamA_acc", "query name": "uniprot_acc", "from":"pfam_seq_start", "to":"pfam_seq_end"})\
                    .assign(pfamA_id=domain_data[pfamA_acc]["pfamA_id"])]
        )

    # Add clan_id
    clan_data = {pfamA_acc:domain_data[pfamA_acc]["clan_acc"] for pfamA_acc in domain_data}
    p_df["clan_acc"] = p_df["pfamA_acc"].map(clan_data)
    clan_accs = p_df["clan_acc"].unique()
    clan_df = clanAccToId(clan_accs)
    # Now add clan_id
    p_df = pd.merge(p_df, clan_df, on="clan_acc", how="outer")
    # Keep only clan_id
    del p_df["clan_acc"]

    # Put in table. Consider only isopeptide bonds with no domain annotation
    t_df = pd.merge(hmm_df[hmm_df["is_domain"]==False][['HMM', 'uniprot_acc', 'hmm_start',
             'hmm_end', 'seq_start', 'seq_end',
            'seq_evalue', 'hit_evalue', 'hit_bitscore', 'Cna', 'taxonomy',
       'species', 'kingdom']].reset_index(), p_df, on=["uniprot_acc"]).set_index("index")
    t_df.index.name = None

    # Assign domain
    def isDomain(row):
        """

            Assign if at least 2/3 residues are in domain

        """
        seq_start = row["seq_start"]
        seq_end = row["seq_end"]
        pfam_seq_start = row["pfam_seq_start"]
        pfam_seq_end = row["pfam_seq_end"]

        matched_domain_len = seq_end-seq_start
        overlap = 0
        for aa in range(seq_start, seq_end):
            if aa >= pfam_seq_start and aa < pfam_seq_end:
                overlap += 1
        status = False

        # Overlap should cover at least 2/3 of the matched domain
        if overlap / matched_domain_len >= 2 / 3:
            status = True
        return status

    t_df["is_domain"] = t_df.apply(isDomain, axis=1)
    t_df = t_df.sort_values("is_domain", ascending=False)\
            .drop_duplicates(["HMM", "uniprot_acc", "seq_start", "seq_end"], keep="first")

    hmm_df.loc[t_df.index] = t_df

    hmm_df.loc[hmm_df["is_domain"]==False, "pfamA_acc"] = np.NaN
    hmm_df.loc[hmm_df["is_domain"]==False, "pfamA_id"] = np.NaN
    hmm_df.loc[hmm_df["is_domain"]==False, "clan_id"] = np.NaN
    hmm_df.loc[hmm_df["is_domain"]==False, "pfam_seq_start"] = np.NaN
    hmm_df.loc[hmm_df["is_domain"]==False, "pfam_seq_end"] = np.NaN

    hmm_df[['HMM', 'uniprot_acc', 'hmm_start', 'hmm_end', 'seq_start', 'seq_end',
       'seq_evalue', 'hit_evalue', 'hit_bitscore', 'Cna', 'taxonomy',
       'species', 'kingdom', 'pfamA_acc', 'pfam_seq_start', 'pfam_seq_end',
       'is_domain', 'pfamA_id', 'clan_id']].to_csv(OUTPUT_TABLE, index=False)
    
def clanAccToId(clan_accs:list) -> pd.DataFrame:
    """

        From clan_acc to clan_id
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    clan_accs = ",".join([f"'{u}'" for u in clan_accs])
    cursor.execute(f"SELECT clan_acc, clan_id \
                     FROM {PFAM_VERSION}.clan \
                     WHERE  clan_acc IN ({clan_accs})")

    output = cursor.fetchall()

    return(
        pd.DataFrame(output, columns=["clan_acc", "clan_id"])
          )

def getSequences(uniprot_accs:list) -> pd.DataFrame:
    """

        Download sequences from pfam
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    uniprot_accs = ",".join([f"'{u}'" for u in uniprot_accs])
    cursor.execute(f"SELECT uniprot_acc, sequence\
                     FROM {PFAM_VERSION}.uniprot \
                     WHERE uniprot_acc IN ({uniprot_accs})")

    output = cursor.fetchall()

    df = pd.DataFrame(output, columns=["uniprot_acc", "sequence"])
    df["sequence"] = df["sequence"].str.decode("utf-8")
    
    return df
    
if __name__ == "__main__":
    main()