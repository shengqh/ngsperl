resultFile<-outFile
mappingResultFile<-parFile1
databaseLogFile<-parFile2
#resultFile<-commandArgs()[7]
#mappingResultFile<-commandArgs()[8]
#databaseLogFile<-commandArgs()[9]

#resultFile<-"3018-KCV-52_53_54.group4Mapping.Result"
#mappingResultFile<-"/scratch/cqs/zhaos/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep/bowtie1_fungus_group4_pm_table/result/fungus_group4_pm_3018-KCV-52_53_54.count"
#databaseLogFile<-"/scratch/cqs/zhaos/vickers/reference/bacteria/group4/20160210.log"

library(ggplot2)
library(reshape)

aggregateCountTable<-function(x,group,method=sum) {
	result<-aggregate(x, list(factor(group)), method)
	row.names(result)<-result[,1]
	result<-result[,-1]
	return(result)
}
groupPie<-function(x,maxCategory=10) {
	nameToCountForFigure<-na.omit(rev(sort(x)))
	if (length(nameToCountForFigure)>maxCategory) {
		nameToCountForFigure<-c(nameToCountForFigure[1:(maxCategory-1)],Other=sum(nameToCountForFigure[-(1:(maxCategory-1))]))
	}
	pie(nameToCountForFigure,col=rainbow(length(nameToCountForFigure)),main=paste0("Mapped Reads: ",as.integer(sum(nameToCountForFigure))))
}


mappingResult<-read.delim(mappingResultFile,header=T,row.names=1)
databaseLog<-read.delim(databaseLogFile,header=T,as.is=T)
#row.names(mappingResult)<-gsub("\\.\\d+","",row.names(mappingResult))

id2Species<-databaseLog$Species
names(id2Species)<-databaseLog$Id
speciesInMappingResult<-id2Species[row.names(mappingResult)]

mappingResult2Species<-aggregateCountTable(mappingResult,speciesInMappingResult)

#str(mappingResult2Species)
write.csv(mappingResult2Species,paste0(resultFile,".toSpecies.csv"))

for (i in 1:ncol(mappingResult2Species)) {
	png(paste0(colnames(mappingResult2Species)[i],".Species.png"),width=2000,height=1500,res=300)
	par(mar=c(2,9,2,9))
	temp<-as.matrix(mappingResult2Species)
	groupPie(temp[,i])
	dev.off()
}

#find the most significant row and combine the result
sigNum<-5
if (nrow(mappingResult2Species)>sigNum) {
	temp<-apply(mappingResult2Species,2,function(x) rev(order(x))[1:sigNum])
	mappingResult2SpeciesSelected<-mappingResult2Species[unique(as.vector(temp)),]
	mappingResult2SpeciesSelected<-rbind(mappingResult2SpeciesSelected,Other=colSums(mappingResult2Species[-unique(as.vector(temp)),,drop=FALSE]))
} else {
	mappingResult2SpeciesSelected<-mappingResult2Species
}

mappingResult2SpeciesSelectedForFigure<-mappingResult2SpeciesSelected
mappingResult2SpeciesSelectedForFigure$Species<-row.names(mappingResult2SpeciesSelected)
mappingResult2SpeciesSelectedForFigure<-melt(mappingResult2SpeciesSelectedForFigure)
colnames(mappingResult2SpeciesSelectedForFigure)<-c("Species","Sample","Reads")

png(paste0(resultFile,".top",sigNum,".png"),width=3000,height=1500,res=300)
ggplot(mappingResult2SpeciesSelectedForFigure,aes(x=Sample,y=Reads,fill=Species))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()




