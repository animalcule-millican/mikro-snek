#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)

# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_file = args[2]

load(input_file)

if(args[3] == 'fwd'){
    errF <- learnErrors(filtFs, nbases = 1e+10, multithread=TRUE, randomize = TRUE) 
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    save(errF, derepFs, sample.names, file = output_file)
}

if(args[3] == 'rev'){
    errR <- learnErrors(filtRs, nbases = 1e+10, multithread=TRUE, randomize = TRUE) 
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    save(errR, derepRs, sample.names, file = output_file)
}