resultFile<-outFile
lengthFileList<-parSampleFile1

lengthFiles<-read.delim(lengthFileList,as.is=T,header=F)

lengthAll<-NULL
for (lengthFile in lengthFiles) {
	lengthOne<-read.delim(lengthFile,as.is=T,header=T)
	lengthOne$Sample<-lengthFile
	lengthAll<-rbind(lengthAll,lengthOne)
}

