library(scRNABatchQC)

filelist1<-read.table("fileList1.txt", header=F, stringsAsFactor=F)
  
result<-scRNABatchQC(inputs=filelist1$V1, 
                     names=filelist1$V2,
                     organism=webgestalt_organism,
                     outputFile=paste0(outFile, ".html"))