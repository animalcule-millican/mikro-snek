#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)
# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
#fwd_input = args[1]
#rev_input = args[2]
inputs = args[1]
output_file = args[2]
input_files <- list.files(input_dir)
for (input in inputs) {
    load(input)
}
#load(fwd_input)
#load(rev_input)

mergers <- mergePairs(dadaFsp, derepFsp, dadaRsp, derepRsp, verbose=TRUE)
seqtab = makeSequenceTable(mergers)
save(mergers, seqtab, file = output_file)