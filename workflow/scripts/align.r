#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)
library(Biostrings)
# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_fasta = args[2]
output_file = args[3]

load(input_file)
seqs = Biostrings::DNAStringSet(dada2::getSequences(seqtab.nochim), use.names = TRUE)
names(seqs) = dada2::getSequences(seqtab.nochim)
mult = DECIPHER::AlignSeqs(seqs, processors = 64)
Biostrings::writeXStringSet(mult, output_align)
save(mult, file = output_file)
