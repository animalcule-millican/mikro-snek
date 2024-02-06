#!/usr/bin/env python3
import argparse
import os
from Bio import Entrez, SeqIO
import time
import concurrent.futures
import textwrap
from random import choices

def get_ids(taxa_target, api_key, email):
    gb_dict = {}
    gb_list = []
    Entrez.email = email
    Entrez.api_key = api_key
    search_term = f"16s ribosomal rna[All Fields] AND {taxa_target}[filter]"
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmode = 'text', retmax = 8000000)
    record = Entrez.read(handle)
    for id in record["IdList"]:
        gb_dict[int(id)] = int(id)
    gb_list = list(gb_dict.values())
    return gb_list


def download_genbank(batch, temp_dir):
    attempt = 0
    success = False
    while attempt < 10 and not success:
        attempt += 1
        try:
            handle = Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text")
            success = True
            print("current batch size accepted")
            wait = int(choices(range(1,4), k=1))
            time.sleep(wait)
        except HTTPError as err:
            if err == 429:
                print("Too many requests, waiting 20 seconds")
                time.sleep(20)
                continue
            else:
                raise
    for rec in SeqIO.parse(handle, "genbank"):
        if len(rec.seq) > 700:
            with open(f"{temp_dir}/{rec.id}.gbk", 'w') as output_handle:
                SeqIO.write(rec, output_handle, 'genbank')

def cat_genbank(temp_dir, output_gbk):
     os.system(f"cat {temp_dir}/*.gbk > {output_gbk}")
     os.system(f"rm {temp_dir}/*.gbk")

def arg_parser():
    parser = argparse.ArgumentParser(description='Download 16s sequences from NCBI')
    parser.add_argument('-b', '--batch', type = int, default=1000, help='Batch size for batched downloads')
    parser.add_argument('-o', '--output', type=str, help='Output directory', required = True)
    parser.add_argument('-t', '--taxa', type=str, help='Which taxa to download files for `bacteria` (default) or `archaea`', required = True, default = "bacteria")
    parser.add_argument('--api', type=str, help='NCBI API key', default = os.environ['NCBI_API_KEY'])
    parser.add_argument('--email', type=str, help='Email address', default = os.environ['NCBI_EMAIL'])
    return parser.parse_args()

def main():
    args = arg_parser()
    temp_dir = f"{args.output}/temp"
    if args.taxa == "bacteria":
        taxa_target = 'bacteria'
        output_file = f"{args.output}/bacterial_16s.gbk"
    elif args.taxa == "archaea":
        taxa_target = 'archaea'
        output_file = f"{args.output}/archaea_16s.gbk"
    else:
        print("Please specify bacteria or archaea")
    id_list = get_ids(taxa_target, args.api, args.email)    
    batch_size = args.batch
    batches = [id_list[i:i+batch_size] for i in range(0, len(id_list), batch_size)]
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(download_genbank, batches, temp_dir)
    
    cat_genbank(temp_dir, output_file)
    if os.path.exists(output_file):
        os.rmdir(temp_dir)
        print(textwrap.dedent(r"""
          ____          Genbank Download Complete.............
          \__/         # ##  
         `(  `^=_ p _###_
          c   /  )  |   /
   _____- //^---~  _c  3
 /  ----^\ /^_\   / --,-
(   |  |  O_| \\_/  ,/
|   |  | / \|  `-- /
(((G   |-----|
      //-----\\
     //       \\
   /   |     |  ^|
   |   |     |   |
   |____|    |____|
  /______)   (_____\
    """)
)

if __name__ == "__main__":
    main()