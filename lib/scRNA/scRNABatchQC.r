rm(list=ls()) 
outFile='PEO4'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''
webgestalt_organism='hsapiens'

setwd('/nobackup/h_vangard_1/shengq2/guoyan/20231004_pipseq_10x_comparison/20231005_scRNA/scRNABatchQC/result')

### Parameter setting end ###

library(scRNABatchQC)

filelist1<-read.table("fileList1.txt", header=F, stringsAsFactor=F)
  
result<-scRNABatchQC(inputs=filelist1$V1, 
                     names=filelist1$V2,
                     organism=webgestalt_organism,
                     createReport=FALSE)

saveRDS(result, paste0(outFile, ".rds"))

generateReport(result$sces,
               result$scesMerge, 
               outputFile=paste0(outFile, ".html"))

#writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
