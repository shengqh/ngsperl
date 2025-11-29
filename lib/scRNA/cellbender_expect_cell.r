rm(list=ls()) 
sample_name='DM_1'
outFile='DM_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20251128_9061_scRNA_mm10_cellbender_redo/cellbender.ratio0.8_01_expect_cells/result/DM_1')

### Parameter setting end ###

source("reportFunctions.R")
source("scRNA_func.r")
flist<-fread(parSampleFile1, header=FALSE)
fpath=flist$V1[1]

cat("Reading", fpath,  "...\n")
obj<-read_scrna_data(fpath)
cat(nrow(obj$counts), "genes and", ncol(obj$counts), "cells.\n")

ratio=as.numeric(read_file_map(parSampleFile2, do_unlist=FALSE)$ratio)
expect_cell=as.integer(ncol(obj$counts) * ratio)
cat("With ratio", ratio, ", final expected cells:", expect_cell, "\n")

writeLines(as.character(expect_cell), paste0(outFile, ".num_cells.txt"))

