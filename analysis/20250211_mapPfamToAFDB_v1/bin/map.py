#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Map Pfam domains known to have Isopeptide bonds to the AFDB
    `conda activate isopeptide`
    `sbatch --job-name="map_pfam_to_afdb" -t 32:00:00 --mem 6GB -e map.err \
            -o map.out --mail-type=END --ntasks 1 --cpus-per-task=48 \
            --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/map.py"`
    Francesco Costa 2025-11-02 fcosta@ebi.ac.uk

"""
import pandas as pd
import glob
import os
from mysql import connector
from pathlib import Path
import sys
sys.path.append("../../bin")
from pfamenv import PFAM_USER, PFAM_HOST, PFAM_PASSWORD, PFAM_PORT, PFAM_VERSION
import re
import numpy as np
from datetime import datetime
from dotenv import load_dotenv
from Bio import SeqIO
import subprocess

load_dotenv("../../.env")

# General
FULL_AFDB_SEQUENCES = os.getenv("FULL_AFDB_SEQUENCES")
# Only domains fetched from Pfam are cached; HMM domains are recomputed
USE_CACHE = True
OUTPUT = "output/afdb_isopep_domains.csv"

# Pfam
AFDB_JESS_SCAN_TABLE = os.getenv("AFDB_JESS_SCAN_TABLE")
PFAM_MAPPED_CACHED = "tmp/pfam_mapped.csv"

# HMM
TMP_HMM_DIR = "tmp/HMM"
HMMSEARCH_OUTPUT = "tmp/hmmsearch.out"
HMMSEARCH = "/hps/software/users/agb/research/francesco/software/hmmer-3.3.2/bin/hmmsearch"
HMMPRESS = "/hps/software/users/agb/research/francesco/software/hmmer-3.3.2/bin/hmmpress"
# Add all additional domains
DOMAINS = ["PF24346", "PF24514", "PF24517", "PF24547", "PF24558", "PF24593", "PF24595", "PF25546", "PF25548", "PF25549", "PF25551", "PF25564"]
HMM_MAPPED_CACHED = "tmp/hmm_df.csv"

def main():
    
    # Get AFDB table
    #afdb_df = pd.DataFrame([i.id for i in SeqIO.parse(FULL_AFDB_SEQUENCES, "fasta")],
    #                        columns=["uniprot_acc"])

    if not os.path.exists(PFAM_MAPPED_CACHED) or not USE_CACHE:
        # For domains not added after pfam_37_0: 
        str_df = pd.read_csv(AFDB_JESS_SCAN_TABLE)
        ibd_domains = str_df[(str_df["probability"]>.65)&(~str_df["taxonomy"].isna())]["pfamA_acc"].unique()
        # Download all sequences matching them
        cnx = connector.connect(user=PFAM_USER,
                                password=PFAM_PASSWORD,
                                port=PFAM_PORT,
                                host=PFAM_HOST)
        cursor = cnx.cursor()
        ibd_domains = ",".join([f"'{u}'" for u in ibd_domains])

        query = f"""SELECT uniprot_acc, uniprot_reg_full.pfamA_acc, pfamA.pfamA_id, seq_start, seq_end \
                    FROM {PFAM_VERSION}.uniprot_reg_full, {PFAM_VERSION}.pfamA \
                    WHERE uniprot_reg_full.pfamA_acc IN ({ibd_domains}) \
                    AND uniprot_reg_full.pfamA_acc=pfamA.pfamA_acc\
                    AND uniprot_reg_full.in_full=1"""

        cursor.execute(query)
        output = cursor.fetchall()
        d_df = pd.DataFrame(output, columns=["uniprot_acc", "pfamA_acc", "pfamA_id", "seq_start", "seq_end"])
        d_df.to_csv("tmp/pfam_domains.csv", index=False)
        
        # Map this to the AFDB
        pfam_mapped_df = pd.merge(afdb_df, d_df, on=["uniprot_acc"])
        pfam_mapped_df.to_csv(PFAM_MAPPED_CACHED, index=False)
    
    pfam_mapped_df = pd.read_csv(PFAM_MAPPED_CACHED)
    # For domains added after pfam_37_0: run hmmscan on the whole AFDB
    hmm_mapped_df = hmm_domain()
    hmm_mapped_df.to_csv(HMM_MAPPED_CACHED, index=False)
    
    # Consider only domains present in pfam
    hmm_mapped_df = pd.merge(hmm_mapped_df, 
        is_in_pfam(hmm_mapped_df["uniprot_acc"].unique()))

    # Give priority to existing domains over hmm obtained ones
    # Exclude hmm domains overlapping with other pfam domains
    over_df = pd.merge(hmm_mapped_df,
            pfam_mapped_df[["uniprot_acc", "pfamA_acc", "seq_start", "seq_end"]].drop_duplicates()\
                .rename(columns={"pfamA_acc":"pfam_pfamA_acc", "seq_start":"pfam_seq_start", "seq_end":"pfam_seq_end"}),
            how="left")

    # Find overlapping domains (more than 30% overlap on the pfam domain)
    over_df = over_df[over_df["overlap"]>.3][["uniprot_acc", "pfamA_acc", "seq_start", "seq_end", "pfamA_id"]].drop_duplicates()
    over_df["overlapping"] = True
    hmm_mapped_df = pd.merge(hmm_mapped_df, over_df, how="left").query('overlapping.isna()')[["uniprot_acc", "pfamA_acc", "seq_start", "seq_end", "pfamA_id"]]

    # Save final output
    pd.concat([pfam_mapped_df, hmm_mapped_df]).drop_duplicates().to_csv(OUTPUT, index=False)

def overlap(row):
    seq_start = row["seq_start"]
    seq_end = row["seq_end"]
    pfam_seq_start = row["pfam_seq_start"]
    pfam_seq_end = row["pfam_seq_end"]

    # Overlap in the two cases
    overlap = 0
    # Case 1: HMM is before than pfam
    if seq_start < pfam_seq_start:
        if seq_end > pfam_seq_start:
            overlap = seq_end - pfam_seq_start
    # Case 2: pfam is before than HMM
    elif pfam_seq_start < seq_start:
        if pfam_seq_end > seq_start:
            overlap = pfam_seq_end - seq_start
    
    # Normalize on target length
    return overlap / (pfam_seq_end - pfam_seq_start) 

def is_in_pfam(uniprot_accs:list) -> pd.DataFrame:
    """
    
        Maps uniprot accs pfam uniprot table
    
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    uniprot_accs = ",".join([f"'{u}'" for u in uniprot_accs])
    cursor.execute(f"SELECT uniprot_acc FROM {PFAM_VERSION}.uniprot \
                WHERE uniprot_acc IN ({uniprot_accs})")

    output = cursor.fetchall()

    return pd.DataFrame(output, columns=["uniprot_acc"])

def hmm_domain() -> pd.DataFrame:
    """
    
        Annotate sequences with new domains created in PFAM using HMMSCAN
    
    """

    # Remove old HMM file
    os.makedirs(TMP_HMM_DIR, exist_ok=True)
    if os.path.isfile(os.path.join(TMP_HMM_DIR, "HMM")):
        os.remove(os.path.join(TMP_HMM_DIR, "HMM"))
    
    domain_data = download_hmm()
    init_hmm()
    hmm_scan()
    p_df = hmm_parse(domain_data)

    return p_df

def hmm_parse(domain_data:dict) -> pd.DataFrame:
    """

    """

    # Parse results
    # NOTE: This code avoids overwriting existing domains <---?
    pf_df = pd.read_csv(HMMSEARCH_OUTPUT, comment='#', 
                          names=["target name", "accession1", "tlen", "query name", "accession2", "qlen", "E-value", "full_seq_score", "full_seq_bias",
                                 "#", "of", "c-Evalue", "i-Evalue", "domain_score", "bias", "hmm_from",
                                 "hmm_to", "ali_from", "ali_to", "from", "to", "acc", "description of target"], sep=r'\s+')
    
    # consider a domain specific GA
    p_df = pd.DataFrame()
    for pfamA_acc in domain_data:
        p_df = pd.concat([p_df, 
                        pf_df[(pf_df["full_seq_score"]>domain_data[pfamA_acc]["sequence_threshold"]) &\
                        (pf_df["domain_score"]>domain_data[pfamA_acc]["domain_threshold"]) &\
                        (pf_df["query name"] == pfamA_acc)]\
                [["target name", "query name", "from", "to"]]\
                    .rename(columns={"target name": "uniprot_acc", "query name": "pfamA_acc", "from":"seq_start", "to":"seq_end"})\
                    .assign(pfamA_id=domain_data[pfamA_acc]["pfamA_id"])]
        )
    
    return p_df

def hmm_scan():
    """

        Run HMM scan

    """
    # Run HMMscan
    cmd = f"{HMMSEARCH} --cpu 48 --domtblout {HMMSEARCH_OUTPUT} {os.path.join(TMP_HMM_DIR, 'HMM')} {FULL_AFDB_SEQUENCES}"
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if r.returncode != 0:
        print("Error encountered during scan")
        print(r.stderr)
        sys.exit()

def init_hmm():
    """

        Init HMM database

    """

    # Init database
    for index_file in glob.glob(os.path.join(TMP_HMM_DIR, "HMM.h*")):
        os.remove(index_file)

    cmd = f"{HMMPRESS} {os.path.join(TMP_HMM_DIR, 'HMM')}"
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if r.returncode != 0:
        print("Error encountered during index creation")
        print(r.stderr)
        sys.exit()

def download_hmm() -> dict:
    """
    
        Download HMM files
    
    """
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
    return domain_data

if __name__ == "__main__":
    main()