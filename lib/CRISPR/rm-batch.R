library(MAGeCKFlute)
cdata <- read.table("rawcount.txt",sep="\t",header=T)
bdata <- read.table("BatchMatrix.txt",sep="\t",header=T)
res <- BatchRemove(mat = cdata, batchMat = bdata)
ccdata <- res$data

write.table(ccdata,file="corrected-count.txt",sep="\t",quote=F,col.names=T,row.names=F)

