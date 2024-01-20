#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)
# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_file = args[2]
orientation = args[3]

load(input_file)

if(orientation == 'fwd'){
    dadaFs <- dada(derepFsp, err=errF, multithread=TRUE, pool = TRUE)
    save(dadaFsp, derepFsp, errF, file = output_file)
}

if(orientation == 'rev'){
    dadaRs <- dada(derepRsp, err=errR, multithread=TRUE, pool = TRUE)
    save(dadaRsp, derepRsp, errR, file = output_file)
}
