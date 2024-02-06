from Bio import SeqIO
import os
import random
from glob import glob
import gzip

def get_seq_metrics(input_dir):
    files = list(glob(os.path.join(input_dir, "*.fast*")))
    # choose a random file from the list
    random_file = random.choice(files)
    full_path = os.path.join(input_dir, random_file)
    file_type = file_format(full_path)
    seq_len, seq_count = seq_count_mean_len(full_path, file_type)
    seqlen = refine_length(seq_len)
    return seq_count, seqlen

def file_format(filename):
    with gzip.open(filename, "rt") as file:
        first_line = file.readline().strip()
        if first_line.startswith(">"):
            return "fasta"
        elif first_line.startswith("@"):
            return "fastq"
    raise ValueError("Unable to determine file format (FASTA or FASTQ)")


def seq_count_mean_len(filename, file_type):
    total_length = 0
    total_count = 0
    with gzip.open(filename, "rt") as file:
        for record in SeqIO.parse(file, file_type):
            total_length += len(record.seq)
            total_count += 1
    mean_len = total_length / total_count if total_count > 0 else 0
    seq_count = total_count // 2 if is_paired_end(filename, file_type) else total_count
    return mean_len, seq_count

def is_paired_end(filename, file_type):
    identifiers = set()
    with gzip.open(filename, "rt") as file:
        for record in SeqIO.parse(file, file_type):
            id_without_end = record.id[:-2]
            if id_without_end in identifiers:
                return True
            identifiers.add(id_without_end)
    return False

def get_chunky(chunks):
    n = chunks
    num_chunks = list(range(1, n+1))
    return num_chunks

def refine_length(seq_len):
    if seq_len < 160:
        result = 150
    elif seq_len < 260:
        result = 250
    elif seq_len > 260:
        result = 300
    return result
