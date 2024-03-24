rm(list=ls()) 
sample_name='CD_Met_01'
outFile='CD_Met_01'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/cellbender_01_expect_cells/result/CD_Met_01')

### Parameter setting end ###

source("reportFunctions.R")
source("scRNA_func.r")

flist<-fread(parSampleFile1, header=FALSE)
fpath=flist$V1[1]

cat("Reading", fpath,  "...\n")
obj<-read_scrna_data(fpath)
cat(nrow(obj$counts), "genes and", ncol(obj$counts), "cells.\n")

writeLines(as.character(ncol(obj$counts)), paste0(outFile, ".num_cells.txt"))
