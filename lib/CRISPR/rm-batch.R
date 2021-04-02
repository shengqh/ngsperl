library(MAGeCKFlute)
cdata <- read.table(parFile1,sep="\t",header=T)
bdata <- read.table(parFile2,sep="\t",header=T)
res <- BatchRemove(mat = cdata, batchMat = bdata)
ccdata <- res$data

write.table(ccdata,file=paste0(outFile, ".corrected-count.txt"),sep="\t",quote=F,col.names=T,row.names=F)

