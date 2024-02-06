#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import os
import random


def file_format(filename):
    with open(filename, "r") as file:
        first_line = file.readline().strip()
        if first_line.startswith(">"):
            return "fasta"
        elif first_line.startswith("@"):
            return "fastq"
    raise ValueError("Unable to determine file format (FASTA or FASTQ)")


def seq_count_mean_len(filename, file_type):
    total_length = 0
    total_count = 0
    for record in SeqIO.parse(filename, file_type):
        total_length += len(record.seq)
        total_count += 1
    mean_len = total_length / total_count if total_count > 0 else 0
    seq_count = total_count // 2 if is_paired_end(filename, file_type) else total_count
    return mean_len, seq_count

def is_paired_end(filename, file_type):
    identifiers = set()
    for record in SeqIO.parse(filename, file_type):
        id_without_end = record.id[:-2]
        if id_without_end in identifiers:
            return True
        identifiers.add(id_without_end)
    return False


def refine_length(seq_len):
    if seq_len < 160:
        result = 150
    elif seq_len < 260:
        result = 250
    elif seq_len > 260:
        result = 300
    return result


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "input", help="Path to directory containing input files")
    parser.add_argument('-c', "config", help="Configuration file for Snakemake workflow")
    args = parser.parse_args()
    return args

def main():
    args = arg_parser()
    files = os.listdir(args.input)
    # choose a random file from the list
    random_file = random.choice(files)
    full_path = os.path.join(args.input, random_file)
    file_type = file_format(full_path)
    seq_len, seq_count = seq_count_mean_len(full_path, file_type)
    seqlen = refine_length(seq_len)
    count_var = f'sequence_count: {seq_count} >> {args.config}'
    len_var = f'sequence_lenght: {seqlen} >> {args.config}'
    os.system(f"echo {count_var}")
    os.system(f"echo {len_var}")

if __name__ == "__main__":
    main()