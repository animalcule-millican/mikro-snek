#!/usr/bin/env Rscript
library(dada2)

# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_tax <- args[2]
output_both <- args[3]
tax_file <- args[4]

load(input_file)
set.seed(8675309) # Initialize random number generator for reproducibility
taxa.species <- addSpecies(taxa, tax_file, tryRC = TRUE, n = 1e10, verbose = TRUE)
taxa.multi <- addSpecies(taxa, tax_file, tryRC = TRUE, n = 1e10, verbose = TRUE, allowMultiple = TRUE)
save(taxa.species, file = output_tax)
save(taxa.species, taxa.multi, file = output_both)