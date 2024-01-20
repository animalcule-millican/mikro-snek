#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)
library(foreach)
library(doParallel)

# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)

# get path to directory with fastq files:
path <- args[1]
# get path to output file
outputF <- args[2]
outputR <- args[3]
# get primer sequences 
FWD = args[4]
REV = args[5]

# get lists of forward and reverse fastq files without extensions and path
Fpat = system(paste0("ls ",path,"/*R1*.fastq.gz"), intern = TRUE)[1] %>% strsplit(., "R1") %>% unlist %>% .[2] %>% paste0("_R1",.)
Rpat = system(paste0("ls ",path,"/*R2*.fastq.gz"), intern = TRUE)[1] %>% strsplit(., "R2") %>% unlist %>% .[2] %>% paste0("_R2",.)

# sort forward and rev files to match
Fs <- sort(list.files(path, pattern=Fpat, full.names = TRUE))
Rs <- sort(list.files(path, pattern=Rpat, full.names = TRUE))

# reverse comp primer sequences
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# create paths to filtered fastq files
Fs.filtN <- file.path(path, "filtN", basename(Fs)) # Put N-filterd files in filtN/ subdirectory
Rs.filtN <- file.path(path, "filtN", basename(Rs))

# filter reads
filterAndTrim(Fs, Fs.filtN, Rs, Rs.filtN, maxN = 0, multithread = TRUE)

# create variable to cutadapt path ### TODO: CHANGE TO A SYSTEM VARIABLE ###
cutadapt <- "/work/adina/millican/.conda/envs/dada2/bin/cutadapt" 

path.cut <- file.path(path, "cutadapt")
# create cutadapt output directory if it doesn't exist
if(!dir.exists(path.cut)) dir.create(path.cut)

# set up cutadapt flags/arguments
Fs.cut <- file.path(path.cut, basename(Fs))
Rs.cut <- file.path(path.cut, basename(Rs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# run cutadapt in parallel
cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)
foreach(i = seq_along(Fs)) %dopar% {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                            "-o", Fs.cut[i], "-p", Rs.cut[i], # output files
                            Fs.filtN[i], Rs.filtN[i])) # input files
}

stopCluster(cl) # stop cluster

# Forward and reverse fastq filenames have the format:
fnFs <- sort(list.files(path.cut, pattern = Fpat, full.names = TRUE))
fnRs <- sort(list.files(path.cut, pattern = Rpat, full.names = TRUE))

# Extract sample names, assuming filenames have format:
sample.names <- unname(sapply(cutFs, get.sample.name))
filtFs = file.path(path.cut, "filtered", basename(cutFs))
filtRs = file.path(path.cut, "filtered", basename(cutRs))

# Filter reads based on quality scores and trim to specified lengths
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                        truncQ = 10, minLen = 50, truncLen=tn.ln, 
                        rm.phix = TRUE, compress = TRUE, multithread = TRUE, matchIDs = TRUE)

save(filtFs, fnFs, sample.names, out, file = outputF)
save(filtRs, fnRs, sample.names, out, file = outputR)