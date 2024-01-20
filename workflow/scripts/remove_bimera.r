#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)
# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file = args[1]
both_seqtabs = args[2]
output_file = args[3]


load(input_file)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)
save(seqtab, seqtab.nochim, file = both_seqtabs)
save(seqtab.nochim, file = output_file)
