#!/usr/bin/env Rscript
library(tidyverse)
library(dada2)

chunk <- function(x, n = 50, count) { 
      # Usage:
      #   chunk(x, n, groups = trunc(len/n))
      # Arguments:
      #   x (vector): 
                  #  The input vector that you want to divide into chunks.
      #   n (integer): 
                  #  The desired number of chunks or groups to split the vector x into.
      #  groups (integer, optional, default = trunc(len/n)): 
                  #  The number of elements in each chunk when evenly dividing the vector x into n chunks. 
                  #  This argument is automatically calculated based on len and n unless explicitly provided.
  groups = trunc(count/n)
  f1 <- as.character(sort(rep(1:n, groups)))
  f <- as.character(c(f1, rep(n, overflow)))
  g <- split(x, f)
  g.names <- names(g)
  g.names.ordered <- as.character(sort(as.numeric(g.names)))
  return(g[g.names.ordered])
}

# import variables passed from snakefile
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
save_dir = paste0(args[2],"/derep_",orientation,)
seq_count <- as.integer(args[3])
n_chunks <- as.integer(args[4]) # optional, defaults to 50
orientation <- args[5]

load(input_file)

if(orientation == 'fwd'){
    result = chunk(derepFs, n_chunks, seq_count)
    
    for (i in names(result)) {
    derepFsp = result[[i]]
        save(derepFsp, errF, file = paste0(save_dir,"_",i,'.RData'))
    }
}

if(orientation == 'rev'){
    result = chunk(derepFs, n_chunks, seq_count)
    
    for (i in names(result)) {
    derepRsp = result[[i]]
        save(derepRsp, errR, file = paste0(save_dir,"_",i,'.RData'))
    }
}

