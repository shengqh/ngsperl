resultFile<-outFile
lengthFileList<-parSampleFile1
groupFileList<-parSampleFile2

library(reshape2)
library(ggplot2)

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")

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

if (groupFileList!="") {
	sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
	#keep the groups with samples in the count table
	sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(lengthAllTable)),]
	
	datBySampleGroup<-mergeTableBySampleGroup(lengthAllTable,sampleToGroup)
	temp<-melt(datBySampleGroup)
	colnames(temp)<-c("Position","Group","Percent")
	png(paste0(resultFile,".png"),width=2000,height=1500,res=300)
	ggplot(temp,aes(x=Position,y=Percent,colour= Group))+geom_line()
	dev.off()
}