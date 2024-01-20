#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pickle
import taxopy

def arg_parser():
    parser = argparse.ArgumentParser(description='Reformat fasta headers to be used with taxonomy assignment in the DADA2 package')
    parser.add_argument('-g', '--genbank', help='Input genbank file', required=False)
    parser.add_argument('-t', '--taxonomy', help='Output fasta file', required=False)
    parser.add_argument('-s', '--species', help='Output fasta file', required=False)
    args = parser.parse_args()
    return args

def build_tax_db():
    taxdb = taxopy.TaxDb(nodes_dmp="/work/adina/millican/.tools/taxdb/nodes.dmp", names_dmp="/work/adina/millican/.tools/taxdb/names.dmp")
    with open('/work/adina/millican/.tools/taxdb/taxdb.pkl', 'wb') as output:
        pickle.dump(taxdb, output, pickle.HIGHEST_PROTOCOL)
    return(taxdb)

def load_tax_db():
    with open("/work/adina/millican/.tools/taxdb/taxdb.pkl", 'rb') as f:
        taxdb = pickle.load(f)
    return taxdb

def parse_tax_ranks(taxa):
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

def get_taxonomy(taxdb, seq_dict):
    for key in seq_dict.keys():
        taxid = seq_dict[key]['taxid']
        ident = seq_dict[key]['accession']
        taxa = taxopy.Taxon(taxid, taxdb)
        tax_dict = {}
        domain, phylum, class_, order, family, genus, species = parse_tax_ranks(taxa)
        tax_dict[ident] = {"accession": ident, "taxid": taxid, "domain": domain, "phylum": phylum, "class": class_, "order": order, "family": family, "genus": genus, "species": species, 'seq': seq_dict[key]['accession']}
    return(tax_dict)

### Now to build the fasta formats needed for DADA2
## for normal tax assignment
#'  >Kingom;Phylum;Class;Order;Family;Genus;   
#'  ACGAATGTGAAGTAA......   

def build_taxonomy_fasta(tax_dict, output_fasta):
    record_list = []
    for key in tax_dict.keys():
        seq = tax_dict[key]['seq']
        name = tax_dict[key]['domain'] + ';' + tax_dict[key]['phylum'] + ';' + tax_dict[key]['class'] + ';' + tax_dict[key]['order'] + ';' + tax_dict[key]['family'] + ';' + tax_dict[key]['genus'] + ';' + tax_dict[key]['species']
        record = SeqRecord(Seq(seq), id=name)
        record_list.append(record)
    with open(output_fasta, 'a') as output_handle:
        SeqIO.write(record_list, output_handle, 'fasta')

## for species assignment
#'  >SeqID genus species  
#'  ACGAATGTGAAGTAA......
def build_species_fasta(tax_dict, output_fasta):
    record_list = []
    for key in tax_dict.keys():
        seq = tax_dict[key]['seq']
        name = f"{tax_dict[key]['accession']} {tax_dict[key]['species']}"
        record = SeqRecord(Seq(seq), id = name)
        record_list.append(record)
    with open(output_fasta, 'a') as output_handle:
        SeqIO.write(record_list, output_handle, 'fasta')

def format_genbank(input_genbank):
    recs = [rec for rec in SeqIO.parse(input_genbank, "genbank")]
    seq_dict = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == 'source']
        for feat in feats:
            taxid = int(feat.qualifiers['db_xref'][0].split(':')[1])
        seq_dict[rec.id] = {'accession': rec.id, 'taxid': taxid, 'seq': rec.seq}
    return(seq_dict)

def main():
    args = arg_parser()
    if os.path.isfile('/work/adina/millican/.tools/taxdb/taxdb.pkl'):
        taxdb = load_tax_db()
    else:
        taxdb = build_tax_db()
    seq_dict = format_genbank(args.genbank)
    tax_dict = get_taxonomy(taxdb, seq_dict)
    build_species_fasta(tax_dict, args.species)
    build_taxonomy_fasta(tax_dict, args.taxonomy)

if __name__ == '__main__':
    main()
