#!/bin/bash
source ~/.bashrc
mamba activate mikro-snek
#source random_directory.sh
DIR_PATH=$1
species=$2
taxonomy=$3

/home/glbrc.org/millican/repos/mikro-snek/workflow/scripts/download_gtdb.py -o $DIR_PATH

if [ -f "$DIR_PATH/refs/gtdb-species.fasta" ] && [ -f "$DIR_PATH/refs/gtdb-taxonomy.fasta" ]; then
    cat $species $DIR_PATH/refs/gtdb-species.fasta > $DIR_PATH/refs/add_species_reference.fna
    cat $taxonomy $DIR_PATH/refs/gtdb-taxonomy.fasta > $DIR_PATH/refs/taxonomy_reference.fna 
    gzip $DIR_PATH/refs/taxonomy_reference.fna 
    gzip $DIR_PATH/refs/add_species_reference.fna
    rm -r $DIR_PATH/gtdb
    rm $DIR_PATH/ref/gtdb-species.fasta $DIR_PATH/ref/gtdb-taxonomy.fasta
else
    echo "$(ls $DIR_PATH/refs/gtdb-species.fasta)"
    echo "$(ls $DIR_PATH/refs/gtdb-taxonomy.fasta)"
fi


