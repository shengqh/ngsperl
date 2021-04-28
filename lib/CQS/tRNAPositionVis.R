options(bitmapType='cairo')

resultFile<-outFile
tRNAPositionFileList<-parSampleFile1
groupFileList<-parSampleFile2
tRNASigFileList<-parSampleFile3
#resultFile<-commandArgs()[7]
#tRNAPositionFileList<-commandArgs()[8]
#groupFileList<-commandArgs()[9]
#tRNASigFileList<-commandArgs()[10]

#setwd("/scratch/cqs/zhaos/temp/testRMoudle/tRNAPositionVis/result")
#resultFile<-"test.result"
#tRNAPositionFileList<-"fileList1.txt"
#groupFileList<-"fileList2.txt"

facetColCount=getFacetColCount(groupFileList)

library(ggplot2)
library(grid)
library(RColorBrewer)
print("Doing tRNA position visualization")
maxPos<-100

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
sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)

#All tRNA position distribution
positionRawAllSamples<-NULL
for (i in 1:nrow(sampleTotRNAPositionFile)) {

	positionRaw<-read.delim(sampleTotRNAPositionFile[i,1],header=T,as.is=T,check.names=F)
	bamInfoFile<-gsub(".tRNA.position$",".info",sampleTotRNAPositionFile[i,1])
	bamInfo<-read.delim(bamInfoFile,header=F,as.is=T,row.names=1,check.names=F)
	totalReads<-as.integer(bamInfo["MappedReads",1]) #normlize by total mapped reads
	
	#positionRaw$Group<-sampleToGroup[row.names(sampleTotRNAPositionFile)[i],1]
	if (! (row.names(sampleTotRNAPositionFile)[i] %in% sampleToGroup[,1])) {
		next;
	}
	groupNames<-sampleToGroup[which(sampleToGroup[,1]==row.names(sampleTotRNAPositionFile)[i]),2]
	positionRaw<-cbind(positionRaw,Group=rep(groupNames,each=nrow(positionRaw)))
	
	positionRaw$Sample<-row.names(sampleTotRNAPositionFile)[i]
	
	positionRaw$absCount<-positionRaw$Count*positionRaw$Percentage
	positionRaw$CountPercentage<-positionRaw$absCount/totalReads  #normlize by total mapped reads
	
	temp<-grep("Undet|Pseudo",positionRaw$Feature)
	if (length(temp)>0) {
		positionRaw<-positionRaw[-temp,]
	}
	
	positionRawAllSamples<-rbind(positionRawAllSamples,positionRaw)
}
positionRawAllSamples$Group<-as.character(positionRawAllSamples$Group)
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
write.csv(positionRawAllSamplestRNAMeanSample,paste0(resultFile,".alltRNAPosition.csv"))
temp<-grep("^[A-Z][a-z][a-z]",positionRawAllSamplestRNAMeanSample$tRNA)
positionRawAllSamplestRNAMeanSample3<-positionRawAllSamplestRNAMeanSample[temp,]
write.csv(positionRawAllSamplestRNAMeanSample3,paste0(resultFile,".all3tRNAPosition.csv"))

#positionRawAllSamplestRNAMeanSample3<-read.csv(paste0(resultFile,".all3tRNAPosition.csv"), row.names = 1)
positionRawAllSamplestRNAMeanSample3<-positionRawAllSamplestRNAMeanSample3[which(positionRawAllSamplestRNAMeanSample3$Position<=maxPos),]
png(paste0(resultFile,".alltRNAPosition.png"),width=3000,height=max(2000, 1000 * length(unique(positionRawAllSamplestRNAMeanSample3$Group))),res=300)
m <- ggplot(positionRawAllSamplestRNAMeanSample3, aes(x = Position,y=CountPercentage,fill=tRNA)) +
  geom_bar(stat="identity")+facet_grid(Group ~ .) +
  theme(legend.key.size = unit(0.4, "cm"))+
  ylab("cumulative read fraction (read counts/total reads)")+
  theme(text = element_text(size=20))+theme(legend.text = element_text(size=16))+
  guides(fill= guide_legend(ncol = 1,keywidth=1, keyheight=1.5))+
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(length(unique(positionRawAllSamplestRNAMeanSample3$tRNA)))) + 
  xlim(-10, maxPos+5)
print(m)
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
print("Doing tRNA position visualization for significant tRNAs")
tRNASigNum<-10
if (is.na(tRNASigFileList) | tRNASigFileList=="") {
	temp<-tapply(positionRawAllSamplesMeanSample$CountPercentage,positionRawAllSamplesMeanSample$Feature,sum)
	tRNASigNames<-names(rev(sort(temp)))[1:tRNASigNum]
	tRNASigFileName<-".highesttRNAPosition.pdf"

	temp<-positionRawAllSamplesMeanSample[which(positionRawAllSamplesMeanSample$Feature %in% tRNASigNames),]
	write.csv(positionRawAllSamplesMeanSample,paste0(resultFile,".tRNAPositionEach.csv"))
	m <- ggplot(temp, aes(x = Position,y=CountPercentage))
	pdf(paste0(resultFile,tRNASigFileName),height=15,width=7)
	print(
			m + geom_bar(stat="identity")+facet_grid(Feature ~ Group)+
			ylab("Read fraction (read counts/total reads)")+
			theme(strip.text.y = element_text(size = 4))
	)
	dev.off()
	
} else {
	tRNASigFiles<-read.delim(tRNASigFileList,as.is=T,header=F,row.names=2)
	tRNASigNames<-NULL
	for (i in 1:nrow(tRNASigFiles)) {
		tRNASigFileEach<-tRNASigFiles[i,1]
		if (file.exists(tRNASigFileEach)) {
			tRNASig<-read.csv(tRNASigFileEach,header=T,row.names=1)
		} else {
			tRNASig<-matrix(NA,ncol=0,nrow=0)
		}
		if (nrow(tRNASig)==0) {
			print(paste0("No significant changed tRNA. Will plot ",tRNASigNum," tRNAs with highest reads"))
			temp<-tapply(positionRawAllSamplesMeanSample$CountPercentage,positionRawAllSamplesMeanSample$Feature,sum)
			tRNASigNames<-names(rev(sort(temp)))[1:tRNASigNum]
			tRNASigFileName<-paste0(".",row.names(tRNASigFiles)[i],".highesttRNAPosition.pdf")
		} else {
			tRNASigNameEach<-row.names(tRNASig)
			tRNASigNameEach<-sapply(strsplit(tRNASigNameEach,";"),function(x) x[1])
			if (length(tRNASigNameEach)>tRNASigNum) {tRNASigNameEach<-tRNASigNameEach[1:tRNASigNum]}
			tRNASigNames<-tRNASigNameEach
			tRNASigFileName<-paste0(".",row.names(tRNASigFiles)[i],".significanttRNAPosition.pdf")
		}
		
		temp<-positionRawAllSamplesMeanSample[which(positionRawAllSamplesMeanSample$Feature %in% tRNASigNames),]
		write.csv(positionRawAllSamplesMeanSample,paste0(resultFile,"tRNAPositionEach.csv"))
		if (nrow(temp)!=0) { #==0 means all significant tRNA has no position information (means they are not high abudance tRNA, not extracted position information by smallRNA count)
			m <- ggplot(temp, aes(x = Position,y=CountPercentage))
			pdf(paste0(resultFile,tRNASigFileName),height=15,width=7)
			print(
					m + geom_bar(stat="identity")+facet_grid(Feature ~ Group)+
							ylab("Read fraction (read counts/total reads)")+
							theme(strip.text.y = element_text(size = 4))
			)
			dev.off()
		}
	}
}

#significant tRNA for each sample
geneselected<-intersect(positionRawAllSamples$Feature,paste0("tRNA:",tRNASigNames))
if (length(geneselected)>0) {
	for (i in 1:length(geneselected)) {
		temp<-positionRawAllSamples[which(positionRawAllSamples$Feature==geneselected[i]),]
		fileSize<-max(length(unique(temp$Sample))/25*2000,2000)
		m <- ggplot(temp, aes(x = Position,ymin=0,ymax=CountPercentage)) + 
				geom_ribbon() +
				ylab("Read fraction (read counts/total reads)") +
				theme(strip.text.y = element_text(size = 4))
		if(facetColCount > 0){
			m <- m + facet_wrap(~Sample, ncol=facetColCount)
		}else{
			m <- m + facet_wrap(~Sample)
		}
		png(paste0(resultFile,tRNASigFileName,".",gsub("tRNA:","",geneselected[i]),".png"),height=fileSize,width=fileSize,res=300)
		print(m)
		dev.off()
	}
}

#Plot position of interested tRNA
if (exists("tRnaPlotPosition")) {
	tRnaPlotPositionMatchData<-NULL
	for (x in tRnaPlotPosition) {
		temp<-unique(positionRawAllSamples$Feature[grep(x,positionRawAllSamples$Feature)])
		if (length(temp)>0) {
			tRnaPlotPositionMatchData<-c(tRnaPlotPositionMatchData,temp)
		}
	}
	tRnaPlotPositionMatchData<-unique(tRnaPlotPositionMatchData)
	if (length(tRnaPlotPositionMatchData)>0) {
		for (i in 1:length(tRnaPlotPositionMatchData)) {
			temp<-positionRawAllSamples[which(positionRawAllSamples$Feature==tRnaPlotPositionMatchData[i]),]
			fileSize<-max(length(unique(temp$Sample))/25*2000,2000)
			m <- ggplot(temp, aes(x = Position,ymin=0,ymax=CountPercentage)) + 
					geom_ribbon() +
					ylab("Read fraction (read counts/total reads)") +
					theme(strip.text.y = element_text(size = 4))
			if(facetColCount > 0){
				m <- m + facet_wrap(~Sample, ncol=facetColCount)
			}else{
				m <- m + facet_wrap(~Sample)
			}
			png(paste0(resultFile,".tRNAInterested.",gsub("tRNA:","",tRnaPlotPositionMatchData[i]),".png"),height=fileSize,width=fileSize,res=300)
			print(m)
			dev.off()
		}
	}
}
