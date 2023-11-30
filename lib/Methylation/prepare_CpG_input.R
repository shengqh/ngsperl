rm(list=ls()) 
sample_name='C003_Control'
outFile='C003_Control'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20231024_10473_WGBS/MethylKitPreparation/result/C003_Control')

### Parameter setting end ###

library(R.utils)
library(data.table)

#require(tidyverse)
dnmtools_path <- paste("..", "..", "dnmtools", "result", sep="/")

#prepare to find all specific format files
cpg.all <- read.table(parSampleFile1, sep="\t")

#read in *.meth files
cpg.infile <- cpg.all$V1[1]
id <- sample_name
#format cpg data frame
DT.cpg <- fread(cpg.infile, header = F, stringsAsFactors = FALSE, data.table=FALSE)
colnames(DT.cpg) <- c("chr", "base", "strand", "Type", "freqC", "coverage")
DT.cpg$chrBase <- paste(DT.cpg$chr, DT.cpg$base, sep = ".")
DT.cpg$freqC <- round(DT.cpg$freqC * 100, 3)
DT.cpg$freqT <- 100 - DT.cpg$freqC
DT.cpg$strand <- "F"
DT.cpg <- DT.cpg[,c("chrBase","chr","base","strand","coverage","freqC","freqT")]
cpg_file=paste(id, "CpG", "txt", sep = ".")
write.table(DT.cpg, file = cpg_file, sep = "\t", quote = F, row.names = F)

gzip(cpg_file, overwrite=TRUE)