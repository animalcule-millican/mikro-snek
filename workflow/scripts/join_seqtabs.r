#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)
# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
path = paste0(args[1],"/rdata")
output_file = args[2]

files = list.files(path, pattern = "merged_dada_*", full.names = TRUE)

for(i in 1:length(files)){
    load(files[i])
    assign(paste0("seqtab.",i), seqtab, envir = .GlobalEnv)
}

tabs = ls(pattern = 'seqtab.')
seqtab = dada2::mergeSequenceTables(tables = tabs)

save(seqtab, file = output_file)