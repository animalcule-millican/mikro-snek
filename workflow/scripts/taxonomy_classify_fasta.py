#!/usr/bin/env python3
import argparse
import os
from urllib.error import HTTPError
from Bio import Entrez, SeqIO
import time
import concurrent.futures
import textwrap
import taxopy
import mmap
from random import choices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle

def build_tax_db():
    taxdb = taxopy.TaxDb(nodes_dmp="/work/adina/millican/.tools/taxdb/nodes.dmp", names_dmp="/work/adina/millican/.tools/taxdb/names.dmp")
    with open('/work/adina/millican/.tools/taxdb/taxdb.pkl', 'wb') as output, open("/work/adina/millican/.tools/taxdb/taxdb.pkl", 'rb') as f:
        pickle.dump(taxdb, output, pickle.HIGHEST_PROTOCOL)
        mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        # Load the pickled object from the memory-mapped file
        taxdb = pickle.load(mmapped_file)
        # Close the memory-mapped file
        mmapped_file.close()
    return(taxdb)

def load_tax_db():
    with open("/work/adina/millican/.tools/taxdb/taxdb.pkl", 'rb') as f:
        # Memory-map the file
        mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        # Load the pickled object from the memory-mapped file
        taxdb = pickle.load(mmapped_file)
        # Close the memory-mapped file
        mmapped_file.close()
    return taxdb

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

def arg_parser():
    parser = argparse.ArgumentParser(description='Download 16s sequences from NCBI')
    parser.add_argument('-b', '--batch', type = int, default=1000, help='Batch size for batched downloads')
    parser.add_argument('-t', '--temp_dir', type=str, help='Temp directory to store files', required = True)
    parser.add_argument('-o', '--output', type=str, help='Output directory', required = True)
    parser.add_argument('--bacteria', action='store_true', help='Download Bacterial 16s sequences from NCBI')
    parser.add_argument('--archaea', action='store_true', help='Download Archaea 16s sequences from NCBI')
    parser.add_argument('--api', type=str, help='NCBI API key', default = "1f5bf9ae3cc3fd00c7fb07b10597be445809")
    parser.add_argument('--email', type=str, help='Email address', default = "millican.m@gmail.com")
    return parser.parse_args()


def parse_tax_ranks(taxid, taxdb):
    taxa = taxopy.Taxon(taxid, taxdb)
    try:
        domain = taxa.rank_name_dictionary['superkingdom']
    except KeyError:
        domain = "NA"
    try:
        phylum = taxa.rank_name_dictionary['phylum']
    except KeyError:
        phylum = "NA"
    try:
        class_ = taxa.rank_name_dictionary['class']
    except KeyError:
        class_ = "NA"
    try:
        order = taxa.rank_name_dictionary['order']
    except KeyError:
        order = "NA"
    try:
        family = taxa.rank_name_dictionary['family']
    except KeyError:
        family = "NA"
    try:
        genus = taxa.rank_name_dictionary['genus']
    except KeyError:
        genus = "NA"
    try:
        species = taxa.rank_name_dictionary['species']
    except KeyError:
        species = "NA"
    return(domain, phylum, class_, order, family, genus, species)

def download_genbank(batch, taxdb):
    attempt = 0
    success = False
    while attempt < 10 and not success:
        attempt += 1
        try:
            handle = Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text")
            success = True
            wait = int(choices(range(1,4))[0])
            time.sleep(wait)
        except HTTPError as err:
            if err == 429:
                print("Too many requests, waiting 20 seconds")
                time.sleep(20)
                continue
            else:
                raise
    recs = [rec for rec in SeqIO.parse(handle, "genbank")]
    seq_dict = {}
    record_list = []
    sp_record = []
    for rec in recs:
        # extract features from source
        feats = [feat for feat in rec.features if feat.type == 'source']
        for feat in feats:
            # parse out taxid from features
            taxid = int(feat.qualifiers['db_xref'][0].split(':')[1])
        
        # get taxonomy lineage information
        #taxa = taxopy.Taxon(taxid, taxdb)
        domain, phylum, class_, order, family, genus, species = parse_tax_ranks(taxid, taxdb)
        
        # Create a seqRecord for the taxonomy fasta
        name = domain + ';' + phylum + ';' + class_ + ';' + order + ';' + family + ';' + genus + ';' + species
        record = SeqRecord(Seq(rec.seq), id=name)
        record_list.append(record)
        
        # If species is assigned, the create a seqRecord for the add species fasta
        if not species is None:
            sp_name = rec.id + ' ' + genus + ' ' + species
            sp_record.append(SeqRecord(Seq(rec.seq), id=sp_name))
        # build a dict for seqRecord information 
        seq_dict[rec.id] = {"accession": rec.id, "taxid": taxid, "domain": domain, "phylum": phylum, "class": class_, "order": order, "family": family, "genus": genus, "species": species, 'seq': rec.seq}
    return(seq_dict, record_list, sp_record)

def write_taxonomy(output_file, record_list):
    with open(output_file, 'a') as output_handle:
        SeqIO.write(record_list, output_handle, 'fasta')

def write_species(sp_output_file, sp_record):
    with open(sp_output_file, 'a') as output_handle:
        SeqIO.write(sp_record, output_handle, 'fasta')

def write_seq_dict(seq_dict_file, combined_seq_dict):
    with open(seq_dict_file, 'wb') as f:
        pickle.dump(combined_seq_dict, f)

def main():
    args = arg_parser()
    if args.bacteria:
        taxa_target = 'bacteria'
        output_file = f"{args.output}/bacteria_taxonomy.fna"
        sp_output_file = f"{args.output}/bacteria_species.fna"
        seq_dict_file = f"{args.output}/bacteria_seq_dict.pkl"
    elif args.archaea:
        taxa_target = 'archaea'
        output_file = f"{args.output}/archaea_taxonomy.fna"
        sp_output_file = f"{args.output}/archaea_species.fna"
        seq_dict_file = f"{args.output}/archaea_seq_dict.pkl"
    else:
        print("Please specify bacteria or archaea")
    
    if os.path.isfile('/work/adina/millican/.tools/taxdb/taxdb.pkl'):
        taxdb = load_tax_db()
    else:
        taxdb = build_tax_db()

    id_list = get_ids(taxa_target, args.api, args.email)    
    batch_size = args.batch
    batches = [id_list[i:i+batch_size] for i in range(0, len(id_list), batch_size)]

    with concurrent.futures.ThreadPoolExecutor(max_workers = 64) as executor:
        futures = [executor.submit(download_genbank, batch, taxdb) for batch in batches]
        # Wait for all batches to complete
        concurrent.futures.wait(futures)
    
    record_list = []
    sp_record = []
    combined_seq_dict = {}
    for future in futures:
        seq_dict, batch_record_list, batch_sp_record = future.result()
        record_list.extend(batch_record_list)
        sp_record.extend(batch_sp_record)
        combined_seq_dict.update(seq_dict)
        
    # submit each function to a separate worker in the thread pool executor
    with concurrent.futures.ThreadPoolExecutor(max_workers = 3) as executor:
        futures = [executor.submit(write_taxonomy, output_file, record_list), 
                   executor.submit(write_species, sp_output_file, sp_record), 
                   executor.submit(write_seq_dict, seq_dict_file, combined_seq_dict)]
    # Wait for all tasks to complete
    concurrent.futures.wait(futures)
    
    if os.path.exists(output_file):
        print(textwrap.dedent(r"""
          ____          Taxonomy ready to go.............
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
