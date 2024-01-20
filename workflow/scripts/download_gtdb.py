#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from Bio import SeqIO # pip install biopython
from Bio.SeqRecord import SeqRecord # pip install biopython
from Bio.Seq import Seq # pip install biopython
from glob import glob # base python package
import os # base python package
import tarfile # base python package
import wget # pip install wget

def download_gtdb(url, output_direcotry):
    if not os.path.exists(output_direcotry):
        os.makedirs(output_direcotry)
    try:
        # Download the file
        filename = wget.download(url, out = output_direcotry)
        if tarfile.is_tarfile(filename):
            with tarfile.open(filename) as f:
                f.extractall(path = output_direcotry)  # Extract all members from the archive to the current working directory
    except Exception as e:
        print(f"An error occurred: {e}")

def fix_taxonomy(taxonomy):
    domain,phylum,_class, order, family, genus, species = taxonomy.split(";")
    domain = domain[3:]
    phylum = phylum[3:]
    _class = _class[3:]
    order = order[3:]
    family = family[3:]
    genus = genus[3:]
    species = species[3:]
    return(domain, phylum, _class, order, family, genus, species)

def read_files(download_path, tax_path, pattern = "*.fna"):
    tax_fasta = f"{tax_path}/gtdb-taxonomy.fasta"
    sp_fasta = f"{tax_path}/gtdb-species.fasta"
    if not os.path.exists(f"{tax_path}"):
        os.makedirs(f"{tax_path}")
    files = list(glob(os.path.join(download_path, pattern)))
    print(files)
    for file in files:
        for record in SeqIO.parse(file, "fasta"):
            seqid = record.id
            header = record.description.split(" ")[1]
            domain, phylum, _class, order, family, genus, species = fix_taxonomy(header)
            record.id = f"{domain};{phylum};{_class};{order};{family};{genus};{species}"
            sp_id = f"{seqid} {genus} {species}"
            seq_record = SeqRecord(Seq(record.seq), id=record.id, description="")
            species_record = SeqRecord(Seq(record.seq), id=sp_id, description="")
            with open(tax_fasta, 'a') as output_handle, open(sp_fasta, 'a') as sp_handle:
                SeqIO.write(seq_record, output_handle, 'fasta')
                SeqIO.write(species_record, sp_handle, 'fasta')

def parse_args():
    parser = argparse.ArgumentParser(description='Download GTDB SSU rRNA sequences')
    parser.add_argument('-o', '--output', help='Path to mikrosnek directory', required=True)
    parser.add_argument('-b', '--bacteria', help='URL to bacterial ssu rep file from GTDB', default="https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/bac120_ssu_reps_r214.tar.gz")
    parser.add_argument('-a', '--archaea', help='URL to bacterial ssu rep file from GTDB', default="https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/ar53_ssu_reps_r214.tar.gz")
    args = parser.parse_args()
    return args

def main():
    # Parse command-line arguments
    args = parse_args()
    download_path = f"{args.output}/gtdb"
    tax_path = f"{args.output}/refs"
    if not os.path.exists(f"{download_path}/bac120_ssu_reps_r214.fna"): 
        # Download GTDB files for bacterial taxa
        download_gtdb(args.bacteria, download_path)
    if not os.path.exists(f"{download_path}/ar53_ssu_reps_r214.fna"): 
        # Download GTDB files for bacterial taxa
        download_gtdb(args.archaea, download_path)
    # process downloaded fasta files in to dada2 formatted files
    if os.path.exists(f"{download_path}/bac120_ssu_reps_r214.fna") and os.path.exists(f"{download_path}/ar53_ssu_reps_r214.fna"):
        read_files(download_path, tax_path) 
    else:
        print("Files did not download correctly")
    if os.path.exists(f"{tax_path}/gtdb-taxonomy.fasta") and os.path.exists(f"{tax_path}/gtdb-species.fasta"):
        print("Taxonomy assignment and add species files created!")
    else:
        print("FAIL!")

if __name__ == '__main__':
    main()