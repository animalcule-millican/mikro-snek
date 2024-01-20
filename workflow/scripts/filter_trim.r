#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)

# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
outputF <- args[2]
outputR <- args[3]
tn.ln <- as.numeric(args[4])

Fpat = system(paste0("ls ",path,"/*R1*.fastq.gz"), intern = TRUE)[1] %>% strsplit(., "R1") %>% unlist %>% .[2] %>% paste0("_R1",.)
Rpat = system(paste0("ls ",path,"/*R2*.fastq.gz"), intern = TRUE)[1] %>% strsplit(., "R2") %>% unlist %>% .[2] %>% paste0("_R2",.)

# sort and list Forward and Reverse files based on patterns determined above
fnFs <- sort(list.files(path, pattern=Fpat, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=Rpat, full.names = TRUE))

# extract samples names - this only uses forward reads.....F and R files should be in same order
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# create a file path and naming for outputting the filtered reads - filterAndTrim function will output here
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter reads based on quality scores and trim to specified lengths
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                        truncQ = 10, minLen = 50, truncLen=tn.ln, 
                        rm.phix = TRUE, compress = TRUE, multithread = TRUE, matchIDs = TRUE)

save(filtFs, fnFs, sample.names, out, file = outputF)
save(filtRs, fnRs, sample.names, out, file = outputR)