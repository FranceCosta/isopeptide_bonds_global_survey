#! /usr/env/python
# -*- coding: utf-8 -*-

""" 

    Get taxonomy, domain annotation, plddt,
    new domain annotations from HMM
    `conda activate isopeptide`
    `sbatch --job-name="annotate_afdb_scan" -t 48:00:00 --mem 16GB -e log/annotate_%j.err \
            -o log/annotate_%j.out --mail-type=END \
            --mail-user=fcosta@ebi.ac.uk --wrap="python3 bin/annotate.py"`
    Francesco Costa 2024-01-08 fcosta@ebi.ac.uk

"""

import pandas as pd
import os
from mysql import connector
from pathlib import Path
import sys
sys.path.append("../../bin")
from pfamenv import PFAM_USER, PFAM_HOST, PFAM_PASSWORD, PFAM_PORT, PFAM_VERSION
import re
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import subprocess
import ast
import numpy as np
from datetime import datetime
import glob 

# General global variables
USE_CACHE = True
CACHE_PLDDT = "tmp/cached_plddt.csv"
CACHE_DOMAIN = "output/ibp_domains.csv"
CACHE_TAX = "tmp/cached_taxonomy.csv"
CACHE_NO_HMM_TABLE = "tmp/cached_no_hmm.csv"
RESULT_GLOB = "/hps/nobackup/agb/research/francesco/tmp/jessAFDB_v7/results/*/jess_out.csv"
DATE = datetime.today().strftime('%Y%m%d')
OUTPUT_TABLE = os.path.join("output", f"{DATE}_adfb_scan.csv")

# HMM domain annotation global variables
CACHE_HMM_DOMAIN = "tmp/cached_hmm_domains.csv"
HMMSCAN = "/hps/software/users/agb/research/francesco/software/hmmer-3.3.2/bin/hmmscan"
HMMPRESS = "/hps/software/users/agb/research/francesco/software/hmmer-3.3.2/bin/hmmpress"
CACHE_SEQS_FASTA = "tmp/cached_sequences.fa"
CACHE_SEQS_TABLE = "tmp/cached_sequences.csv"
TMP_HMM_DIR = "tmp/HMMs"
HMMSCAN_OUTPUT = "tmp/hmscan.out"
CACHE_HMM_TABLE = "tmp/cached_hmm.csv"
DOMAINS = ["PF24346", "PF24514", "PF24517", "PF24547", "PF24558", "PF24593", "PF24595", "PF25546", "PF25548", "PF25549", "PF25551", "PF25564"]

def main():
    
    columns = ("struct_file", "protein_name", "probability", "chain", "r1_bond", 
            "r_cat", "r2_bond", "r1_bond_name", "r_cat_name", "r2_bond_name", 
            "bond_type", "rmsd", "r_asa", "template")

    df = pd.DataFrame(columns=columns)

    for file in glob.glob(RESULT_GLOB):
        tmp_df = pd.read_csv(file, names=columns)
        df = pd.concat([df, tmp_df])
    
    # Remove unnecessary entries
    df = df[df["rmsd"]<=1]

    # Get domains
    print("Running domains")
    uniprot_accs = df["protein_name"].unique()
    if os.path.isfile(CACHE_DOMAIN) == False or not USE_CACHE:
        domain_df = domain(uniprot_accs)
        domain_df.to_csv(CACHE_DOMAIN, index=False)
    else:
        domain_df = pd.read_csv(CACHE_DOMAIN)

    # Get taxonomy
    print("Running taxonomy")
    uniprot_accs = df["protein_name"].unique()
    if os.path.isfile(CACHE_TAX) == False or not USE_CACHE:
        tax_df = taxonomy(uniprot_accs)
        tax_df.to_csv(CACHE_TAX, index=False)
    else:
        tax_df = pd.read_csv(CACHE_TAX)

    # Format table
    print("Formatting table")
    df = pd.merge(df.rename(columns={"protein_name":"uniprot_acc"}), domain_df, on="uniprot_acc", how="left")
    df = pd.merge(df, tax_df, how="left")
    df["kingdom"] = df["taxonomy"].str.replace(".", "").fillna("ND").apply(lambda x: x.split(";")[0]).replace("ND", np.NaN)
    # Delete df with is_domain false IF DUPLICATED
    df["is_domain"] = df.apply(is_domain, axis=1)
    df = df.sort_values("is_domain", ascending=False) \
        .drop_duplicates(["uniprot_acc", "r1_bond", "r_cat", "r2_bond"], keep="first")
    # Save for debugging purposes
    df.to_csv(CACHE_NO_HMM_TABLE, index=False)

    # Get HMM domains and format the table
    print("Running HMM domains")
    df = hmm_domain(df)
    # Clear from unnecessary data
    df.loc[df["is_domain"]==False, "pfamA_acc"] = np.NaN
    df.loc[df["is_domain"]==False, "pfamA_id"] = np.NaN
    df.loc[df["is_domain"]==False, "clan_id"] = np.NaN
    df.loc[df["is_domain"]==False, "clan_acc"] = np.NaN
    df.loc[df["is_domain"]==False, "seq_start"] = np.NaN
    df.loc[df["is_domain"]==False, "seq_end"] = np.NaN
    df.to_csv(CACHE_HMM_TABLE, index=False)

    # Get plddt
    print("Running plddt")
    if os.path.isfile(CACHE_PLDDT) == False or not USE_CACHE:
        plddt_df = plddt(df)
        plddt_df.to_csv(CACHE_PLDDT, index=False)
    else:
        plddt_df = pd.read_csv(CACHE_PLDDT)

    df = pd.merge(df, plddt_df, how="left")
    split = pd.DataFrame(df['plddt'].apply(eval).to_list(), columns = ["r1_plddt", "r2_plddt", "r3_plddt"])
    df = pd.concat([df, split], axis=1)
    del df["plddt"]

    # Save
    df.to_csv(OUTPUT_TABLE, index=False)

def hmm_domain(df:pd.DataFrame) -> pd.DataFrame:
    """
    
        Annotate sequences with new domains created in PFAM using HMMSCAN
    
    """
    df = df.copy()
    # Remove old HMM file

    os.makedirs(TMP_HMM_DIR, exist_ok=True)
    if os.path.isfile(os.path.join(TMP_HMM_DIR, "HMM")):
        os.remove(os.path.join(TMP_HMM_DIR, "HMM"))
    
    # Get sequences that hmmscan will use
    if not os.path.isfile(CACHE_SEQS_FASTA) or not USE_CACHE:
        seq_df = get_sequences(df["uniprot_acc"].unique())
        seq_df.to_csv(CACHE_SEQS_TABLE, index=False)
        with open(CACHE_SEQS_FASTA, "w") as fh:
            for _, row in seq_df.iterrows():
                fh.write(">"+row["uniprot_acc"]+"\n"+row["sequence"]+"\n")
    domain_data = download_hmm()
    init_hmm()
    hmm_scan()
    p_df = hmm_parse(domain_data)
    
    # Add clan_id
    clan_data = {pfamA_acc:domain_data[pfamA_acc]["clan_acc"] for pfamA_acc in domain_data}
    p_df["clan_acc"] = p_df["pfamA_acc"].map(clan_data)
    clan_accs = p_df["clan_acc"].unique()
    clan_df = clanacc_to_id(clan_accs)
    # Now add clan_id
    p_df = pd.merge(p_df, clan_df, on="clan_acc", how="outer")

    # Put in table
    # It is important to keep the original indexing here
    t_df = pd.merge(df[df["is_domain"]==False].reset_index()[[col for col in df.columns if col not in p_df]+["uniprot_acc", "index"]],
                p_df, on="uniprot_acc").set_index("index")
    t_df.index.name = None

    # Assign domain
    t_df["is_domain"] = t_df.apply(is_domain, axis=1)
    t_df = t_df.sort_values("is_domain", ascending=False)\
               .drop_duplicates(["uniprot_acc", "r1_bond", "r_cat", "r2_bond"], keep="first")
    # New data are put back into the table using index
    df.loc[t_df.index] = t_df
    return df

def is_domain(row):
    """

        Assign if at least 2/3 residues are in domain
    
    """
    residues = [row["r1_bond"], row["r_cat"], row["r2_bond"]]
    status = False
    c = 0
    for res in residues:
        if res >= row["seq_start"] and res <= row["seq_end"]:
            c += 1
            #status = True
            #break
    if c >= 2:
        status = True
    return status

def clanacc_to_id(clan_accs:list) -> pd.DataFrame:
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

def hmm_parse(domain_data:dict) -> pd.DataFrame:
    """

    """

    # Parse results
    # NOTE: This code avoids overwriting existing domains (<--??)
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
                   .rename(columns={"target name": "pfamA_acc", "query name": "uniprot_acc", "from":"seq_start", "to":"seq_end"})\
                    .assign(pfamA_id=domain_data[pfamA_acc]["pfamA_id"])]
        )
    
    return p_df

def hmm_scan():
    """

        Run HMM scan

    """
    # Run HMMscan
    cmd = f"{HMMSCAN} --domtblout {HMMSCAN_OUTPUT} {os.path.join(TMP_HMM_DIR, 'HMM')} {CACHE_SEQS_FASTA}"
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
    
def get_sequences(uniprot_accs:list) -> pd.DataFrame:
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

def plddt(df:pd.DataFrame) -> pd.DataFrame:
    """
    
        Extract plDDT once for every protein
    
    """
    df = df[["struct_file", "uniprot_acc", "r1_bond", "r_cat", "r2_bond"]].sort_values(["uniprot_acc", "r1_bond"]).copy()
    outlist = []
    protein_name = None
    for index, row in df.iterrows():
        # Load the new structure only when necessary
        tmp_protein_name = row["uniprot_acc"]
        if tmp_protein_name != protein_name:
            protein_name = tmp_protein_name
            struct_file = row["struct_file"]
            pdb_file = pdb.PDBFile.read(row["struct_file"])
            plddt = pdb_file.get_b_factor()[0]
            structure = pdb_file.get_structure()[0]
            structure.set_annotation("plddt", plddt)
        r1 = row["r1_bond"]
        r2 = row["r_cat"]
        r3 = row["r2_bond"]
        plddt = [atom.plddt for atom in structure if atom.res_id in [r1, r2, r3] and atom.atom_name == "CA"]
        outlist.append(plddt)
    df["plddt"] = outlist

    return df        
    
def domain(uniprot_accs:list) -> pd.DataFrame:
    """
    
        Maps uniprot accs to pfam domains
    
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    uniprot_accs = ",".join([f"'{u}'" for u in uniprot_accs])

    # This query ensures that the domains are added even if a clan annotation is missing!
    query = f"""SELECT 
                    uniprot_acc, 
                    uniprot_reg_full.pfamA_acc, 
                    pfamA.pfamA_id, 
                    clan_membership.clan_acc, 
                    clan.clan_id, 
                    seq_start, 
                    seq_end
                FROM 
                    {PFAM_VERSION}.uniprot_reg_full
                LEFT JOIN 
                    {PFAM_VERSION}.clan_membership 
                    ON uniprot_reg_full.pfamA_acc = clan_membership.pfamA_acc
                LEFT JOIN 
                    {PFAM_VERSION}.clan 
                    ON clan_membership.clan_acc = clan.clan_acc
                JOIN 
                    {PFAM_VERSION}.pfamA 
                    ON uniprot_reg_full.pfamA_acc = pfamA.pfamA_acc
                WHERE 
                    uniprot_acc IN ({uniprot_accs})
                    AND uniprot_reg_full.in_full = 1;"""

    cursor.execute(query)
    output = cursor.fetchall()

    return(pd.DataFrame(output, columns=["uniprot_acc", "pfamA_acc", "pfamA_id", "clan_acc", "clan_id", "seq_start", "seq_end"]))

def taxonomy(uniprot_accs:list) -> pd.DataFrame:
    """
    
        Maps uniprot accs to taxonomy
    
    
    """
    cnx = connector.connect(user=PFAM_USER,
                            password=PFAM_PASSWORD,
                            port=PFAM_PORT,
                            host=PFAM_HOST)
    cursor = cnx.cursor()
    uniprot_accs = ",".join([f"'{u}'" for u in uniprot_accs])
    cursor.execute(f"SELECT uniprot_acc, taxonomy, species FROM {PFAM_VERSION}.uniprot \
                WHERE uniprot_acc IN ({uniprot_accs})")

    output = cursor.fetchall()

    return pd.DataFrame(output, columns=["uniprot_acc", "taxonomy", "species"])

if __name__ == "__main__":
    main()
