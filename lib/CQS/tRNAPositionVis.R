resultFile<-commandArgs()[7]
tRNAPositionFileList<-commandArgs()[8]
groupFileList<-commandArgs()[9]
tRNASigFileList<-commandArgs()[10]

#setwd("/scratch/cqs/zhaos/temp/testRMoudle/tRNAPositionVis/result")
#resultFile<-"test.result"
#tRNAPositionFileList<-"fileList1.txt"
#groupFileList<-"fileList2.txt"

library(ggplot2)
library(grid)

summaryPositionInSamples<-function(positionData,positionKeyVar="PositionKey",
		keepVar=c("Position","Group","tRNA","Count","CountPercentage"),
		summaryVar=c("Count","CountPercentage"),
		summaryFun=mean) {
	positionData<-split(positionData,positionData[,positionKeyVar])
	positionInSamples<-t(sapply(positionData,function(x) {
						result<-x[1,unique(c(keepVar,summaryVar))]
						for (summaryVarEach in summaryVar) {
							result[summaryVarEach]<-summaryFun(x[,summaryVarEach])
						}
						return(result)
					}))
	positionInSamples<-as.data.frame(positionInSamples)
	positionInSamples[]<-lapply(positionInSamples,function(y) if(class(y[[1]])=="numeric"|class(y[[1]])=="integer") as.numeric(y) else as.character(y) )
	return(positionInSamples)
}

sampleTotRNAPositionFile<-read.delim(tRNAPositionFileList,as.is=T,header=F,row.names=2)
sampleToGroup<-read.delim(groupFileList,as.is=T,header=F,row.names=1)

#All tRNA position distribution
positionRawAllSamples<-NULL
for (i in 1:nrow(sampleTotRNAPositionFile)) {

	positionRaw<-read.delim(sampleTotRNAPositionFile[i,1],header=T,as.is=T)
	
	positionRaw$Group<-sampleToGroup[row.names(sampleTotRNAPositionFile)[i],1]
	positionRaw$Sample<-row.names(sampleTotRNAPositionFile)[i]
	
	positionRaw$absCount<-positionRaw$Count*positionRaw$Percentage
	positionRaw$CountPercentage<-positionRaw$absCount/sum(positionRaw$absCount)
	
	temp<-grep("Undet|Pseudo",positionRaw$Feature)
	if (length(temp)>0) {
		positionRaw<-positionRaw[-temp,]
	}
	
	positionRawAllSamples<-rbind(positionRawAllSamples,positionRaw)
}
positionRawAllSamples$tRNA<-substr(sapply(strsplit(positionRawAllSamples$Feature,"-"),function(x) x[2]),0,3)
positionRawAllSamples$PositionKey1<-paste0(positionRawAllSamples$Feature,"_",positionRawAllSamples$Position)
positionRawAllSamples$PositionKey2<-paste0(positionRawAllSamples$tRNA,"_",positionRawAllSamples$Position,"_",positionRawAllSamples$Sample)
positionRawAllSamples$PositionKey3<-paste0(positionRawAllSamples$tRNA,"_",positionRawAllSamples$Position)

positionRawAllSamplestRNAMeanSample<-NULL
for (groupNameEach in unique(positionRawAllSamples$Group)) {
	temp1<-positionRawAllSamples[which(positionRawAllSamples$Group==groupNameEach),]
	temp1<-summaryPositionInSamples(positionData=temp1,positionKeyVar="PositionKey2",
			keepVar=c("Position","Group","Sample","tRNA","PositionKey3"),
			summaryVar=c("CountPercentage"),
			summaryFun=sum)
	temp1<-summaryPositionInSamples(positionData=temp1,positionKeyVar="PositionKey3",
			keepVar=c("Position","Group","tRNA"),
			summaryVar=c("CountPercentage"),
			summaryFun=mean)
	positionRawAllSamplestRNAMeanSample<-rbind(positionRawAllSamplestRNAMeanSample,temp1)
}

m <- ggplot(positionRawAllSamplestRNAMeanSample, aes(x = Position,y=CountPercentage,fill=tRNA))
pdf(paste0(resultFile,".alltRNAPosition.pdf"),width=6,height=6)
m + geom_bar(stat="identity")+facet_grid(Group ~ .) +
		theme(legend.key.size = unit(0.4, "cm"))+
		ylab("cumulative read fraction (read counts/total reads)")
dev.off()

#significant tRNA position distribution
positionRawAllSamplesMeanSample<-NULL
for (groupNameEach in unique(positionRawAllSamples$Group)) {
	temp1<-positionRawAllSamples[which(positionRawAllSamples$Group==groupNameEach),]
	temp1<-summaryPositionInSamples(positionData=temp1,positionKeyVar="PositionKey1",
			keepVar=c("Feature","Position","Group"),
			summaryVar=c("CountPercentage"),
			summaryFun=mean)
	positionRawAllSamplesMeanSample<-rbind(positionRawAllSamplesMeanSample,temp1)
}
positionRawAllSamplesMeanSample$Feature<-gsub("tRNA:","",positionRawAllSamplesMeanSample$Feature)

#significant tRNA names
tRNASigNum<-10
if (is.null(tRNASigFileList)) {
	tRNASigNames<-unique(positionRawAllSamplesMeanSample$Feature)[1:tRNASigNum]
} else {
	tRNASigFiles<-read.delim(tRNASigFileList,as.is=T,header=F,row.names=2)
	tRNASigNames<-NULL
	for (tRNASigFileEach in tRNASigFiles[,1]) {
		tRNASig<-read.csv(tRNASigFileEach,header=T,row.names=1)
		tRNASigNameEach<-row.names(tRNASig)
		tRNASigNameEach<-sapply(strsplit(tRNASigNameEach,";"),function(x) x[1])
		if (length(tRNASigNameEach)>as.integer(tRNASigNum/nrow(tRNASigFiles))) {tRNASigNameEach<-tRNASigNameEach[1:as.integer(tRNASigNum/nrow(tRNASigFiles))]}
		tRNASigNames<-c(tRNASigNames,tRNASigNameEach)
	}
}

temp<-positionRawAllSamplesMeanSample[which(positionRawAllSamplesMeanSample$Feature %in% tRNASigNames),]
m <- ggplot(temp, aes(x = Position,y=CountPercentage))
pdf(paste0(resultFile,".significanttRNAPosition.pdf"),height=15,width=7)
m + geom_bar(stat="identity")+facet_grid(Feature ~ Group)+
		ylab("Read fraction (read counts/total reads)")+
		theme(strip.text.y = element_text(size = 4))
dev.off()


