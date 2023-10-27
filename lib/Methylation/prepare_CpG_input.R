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

#require(tidyverse)
dnmtools_path <- paste("..", "..", "dnmtools", "result", sep="/")

#prepare to find all specific format files
cpg.all <- read.table(parSampleFile1, sep="\t")

#count average bvalue for each hmr of each sample
#read in *.meth files
cpg.infile <- cpg.all$V1[1]
id <- sample_name
#format cpg data frame
DT.cpg <- read.table(cpg.infile, header = F, stringsAsFactors = FALSE)
colnames(DT.cpg) <- c("chr", "base", "strand", "Type", "freqC", "coverage")
DT.cpg$chrBase <- paste(DT.cpg$chr, DT.cpg$base, sep = ".")
DT.cpg$freqC <- round(DT.cpg$freqC * 100, 3)
DT.cpg$freqT <- 100 - DT.cpg$freqC
DT.cpg$strand <- "F"
DT.cpg <- DT.cpg[,c("chrBase","chr","base","strand","coverage","freqC","freqT")]
write.table(DT.cpg, file = paste(id, "CpG", "txt", sep = "."), sep = "\t", quote = F, row.names = F)
