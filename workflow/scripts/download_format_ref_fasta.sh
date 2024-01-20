#!/bin/bash
source ~/.bashrc
export PATH=/work/adina/millican/.tools/mikro-snek/workflow/scripts:$PATH

download_genbank.py -b 3000 -t $TMPDIR -o /work/adina/millican/.tools/mikro-snek --bacteria 
genbank_formatter.py -g /work/adina/millican/.tools/mikro-snek/bacteria_16s.gbk -t /work/adina/millican/.tools/mikro-snek/bacteria_taxonomy.fna -s /work/adina/millican/.tools/mikro-snek/bacteria_species.fna

download_genbank.py -b 3000 -t $TMPDIR -o /work/adina/millican/.tools/mikro-snek --archaea 
genbank_formatter.py -g /work/adina/millican/.tools/mikro-snek/archaea_16s.gbk -t /work/adina/millican/.tools/mikro-snek/archaea_taxonomy.fna -s /work/adina/millican/.tools/mikro-snek/archaea_species.fna



cat /work/adina/millican/.tools/mikro-snek/*_taxonomy.fna > /work/adina/millican/.tools/mikro-snek/NCBI_16s_taxonomy.fna
cat /work/adina/millican/.tools/mikro-snek/*_species.fna > /work/adina/millican/.tools/mikro-snek/NCBI_16s_species_assign.fna

gzip /work/adina/millican/.tools/mikro-snek/NCBI_16s_taxonomy.fna
gzip /work/adina/millican/.tools/mikro-snek/NCBI_16s_species_assign.fna

rm /work/adina/millican/.tools/mikro-snek/bacteria_16s.gbk /work/adina/millican/.tools/mikro-snek/archaea_16s.gbk /work/adina/millican/.tools/mikro-snek/bacteria_taxonomy.fna /work/adina/millican/.tools/mikro-snek/bacteria_species.fna /work/adina/millican/.tools/mikro-snek/archaea_taxonomy.fna /work/adina/millican/.tools/mikro-snek/archaea_species.fna