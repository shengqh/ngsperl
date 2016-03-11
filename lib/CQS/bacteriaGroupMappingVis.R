resultFile<-outFile
mappingResultFile<-parFile1
databaseLogFile<-parFile2
groupFileList<-parSampleFile1
#resultFile<-commandArgs()[7]
#mappingResultFile<-commandArgs()[8]
#databaseLogFile<-commandArgs()[9]

#resultFile<-"3018-KCV-52_53_54.group4Mapping.Result"
#mappingResultFile<-"/scratch/cqs/zhaos/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep/bowtie1_fungus_group4_pm_table/result/fungus_group4_pm_3018-KCV-52_53_54.count"
#databaseLogFile<-"/scratch/cqs/zhaos/vickers/reference/bacteria/group4/20160210.log"

library(ggplot2)
library(reshape)


expandCountTableByName<-function(x,sep=";") {
	namesEach<-strsplit(row.names(x),sep)
	namesEachUnlist<-unlist(namesEach)
	namesEachLength<-sapply(namesEach,length)
	xEach<-x/namesEachLength
	result<-xEach[rep(seq.int(1,nrow(xEach)), namesEachLength),]
	row.names(result)<-namesEachUnlist
	return(result)
}
aggregateCountTable<-function(x,group,method=sum) {
	result<-aggregate(x, list(factor(group)), method)
	row.names(result)<-result[,1]
	result<-result[,-1]
	return(result)
}
groupPie<-function(x,maxCategory=10,main="") {
	nameToCountForFigure<-na.omit(rev(sort(x)))
	if (length(nameToCountForFigure)>maxCategory) {
		nameToCountForFigure<-c(nameToCountForFigure[1:(maxCategory-1)],Other=sum(nameToCountForFigure[-(1:(maxCategory-1))]))
	}
	pie(nameToCountForFigure,col=rainbow(length(nameToCountForFigure)),main=main)
}

mergeTableBySampleGroup<-function(x,sampleToGroup) {
	xRatio<-t(t(x)/colSums(x))
	groupLength<-length(unique(sampleToGroup[,2]))
	xRatioGroupMean<-matrix(NA,ncol=groupLength,nrow=nrow(x))
	colnames(xRatioGroupMean)<-unique(sampleToGroup[,2])
	row.names(xRatioGroupMean)<-row.names(x)
	for (i in 1:groupLength) {
		currentSample<-sampleToGroup[which(sampleToGroup[,2]==colnames(xRatioGroupMean)[i]),1]
		xRatioGroupMean[,i]<-rowMeans(xRatio[,currentSample])
	}
	return(xRatioGroupMean)
}

groupBarplot<-function(x,maxCategory=5,groupName="Species") {
	if (nrow(x)>maxCategory) {
		temp<-apply(x,2,function(y) rev(order(y))[1:maxCategory])
		mappingResult2SpeciesSelected<-x[unique(as.vector(temp)),]
		mappingResult2SpeciesSelected<-rbind(mappingResult2SpeciesSelected,Other=colSums(x[-unique(as.vector(temp)),,drop=FALSE]))
	} else {
		mappingResult2SpeciesSelected<-x
	}
	
	mappingResult2SpeciesSelectedForFigure<-mappingResult2SpeciesSelected
	mappingResult2SpeciesSelectedForFigure$Groups<-row.names(mappingResult2SpeciesSelected)
	mappingResult2SpeciesSelectedForFigure<-melt(mappingResult2SpeciesSelectedForFigure)
	colnames(mappingResult2SpeciesSelectedForFigure)<-c("Groups","Sample","Reads")
	
	ggplot(mappingResult2SpeciesSelectedForFigure,aes(x=Sample,y=Reads,fill=Groups))+
			geom_bar(stat="identity")+
			guides(fill= guide_legend(title = groupName))+
			theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)
databaseLog<-read.delim(databaseLogFile,header=T,as.is=T)
#row.names(mappingResult)<-gsub("\\.\\d+","",row.names(mappingResult))
if (groupFileList!="") {
	sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
	#keep the groups with samples in the count table
	sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(mappingResult)),]
}

id2Species<-databaseLog$Species
names(id2Species)<-databaseLog$Id

mappingResultExpand<-expandCountTableByName(mappingResult)
speciesInMappingResult<-id2Species[row.names(mappingResultExpand)]

mappingResult2Species<-aggregateCountTable(mappingResultExpand,speciesInMappingResult)

#str(mappingResult2Species)
write.csv(mappingResult2Species,paste0(resultFile,".toSpecies.csv"))

for (i in 1:ncol(mappingResult2Species)) {
	png(paste0(resultFile,"_",colnames(mappingResult2Species)[i],".Species.png"),width=2000,height=1500,res=300)
	par(mar=c(2,9,2,9))
	temp<-as.matrix(mappingResult2Species)
	groupPie(temp[,i])
	dev.off()
}

if (groupFileList!="") {
	mappingResult2SpeciesBySampleGroup<-mergeTableBySampleGroup(mappingResult2Species,sampleToGroup)
	
	groupNames<-colnames(mappingResult2SpeciesBySampleGroup)
	for (i in 1:length(groupNames)) {
		png(paste0(resultFile,"_",groupNames[i],".Group.Species.png"),width=2000,height=1500,res=300)
		par(mar=c(2,9,2,9))
		temp<-as.matrix(mappingResult2SpeciesBySampleGroup)
		groupPie(temp[,i],main=paste0("Group: ",groupNames[i]))
		dev.off()
	}
}

#find the most significant row and combine the result
png(paste0(resultFile,".top.png"),width=3000,height=1500,res=300)
groupBarplot(mappingResult2Species,groupName="Species")
dev.off()




