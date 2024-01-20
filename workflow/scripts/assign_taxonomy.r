#!/usr/bin/env Rscript
library(dada2)

# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_tax <- args[2]
output_both <- args[3]
tax_file <- args[4]

set.seed(8675309) # Initialize random number generator for reproducibility
taxa <- assignTaxonomy(seqtab.nochim, tax_file, multithread=TRUE, tryRC = TRUE, taxLevels = c('domain','phylum','class','order','family','genus','species'), verbose=TRUE)
save(taxa, file = output_tax)
save(taxa, seqtab.nochim, file = output_both)