resultFile<-outFile
lengthFileList<-parSampleFile1

library(reshape2)

lengthFiles<-read.delim(lengthFileList,as.is=T,header=F)

lengthAll<-NULL
for (i in 1:nrow(lengthFiles)) {
	lengthFile<-lengthFiles[i,1]
	sampleName<-lengthFiles[i,2]
	lengthOne<-read.delim(lengthFile,as.is=T,header=T)
	lengthOne$Sample<-sampleName
	lengthAll<-rbind(lengthAll,lengthOne)
}

lengthAllTable<-acast(lengthAll,Len~Sample,value.var="Count")
write.csv(lengthAllTable,paste0(resultFile,".csv"))
