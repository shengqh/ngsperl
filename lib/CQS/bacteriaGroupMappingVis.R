resultFile<-commandArgs()[7]
mappingResultFile<-commandArgs()[8]
databaseLogFile<-commandArgs()[9]

#resultFile<-"3018-KCV-52_53_54.group4Mapping.Result"
#mappingResultFile<-"/scratch/cqs/zhaos/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep/bowtie1_fungus_group4_pm_table/result/fungus_group4_pm_3018-KCV-52_53_54.count"
#databaseLogFile<-"/scratch/cqs/zhaos/vickers/reference/bacteria/group4/20160210.log"

library(ggplot2)
library(reshape)

mappingResult<-read.delim(mappingResultFile,header=T,row.names=1)
databaseLog<-read.delim(databaseLogFile,header=T,as.is=T)
row.names(mappingResult)<-gsub("\\.\\d+","",row.names(mappingResult))

id2Species<-databaseLog$Species
names(id2Species)<-databaseLog$Id

mappingResult2Species<-matrix(0,ncol=ncol(mappingResult),nrow=length(unique(id2Species)))
mappingResult2Species<-as.data.frame(mappingResult2Species)
row.names(mappingResult2Species)<-unique(id2Species)
colnames(mappingResult2Species)<-colnames(mappingResult)

for (i in 1:nrow(mappingResult)) {
	temp<-strsplit(row.names(mappingResult)[i],";")[[1]]
	for (z in temp) {
		mappingResult2Species[id2Species[z],]<-as.vector(mappingResult2Species[id2Species[z],])+as.vector(mappingResult[i,]/length(temp))
	}
}

#str(mappingResult2Species)
write.csv(mappingResult2Species,paste0(resultFile,".toSpecies.csv"))

for (i in 1:ncol(mappingResult2Species)) {
	png(paste0(colnames(mappingResult2Species)[i],".bacteriaGroup.png"),width=3000,height=1500,res=300)
	par(mar=c(2,9,2,9))
	temp<-mappingResult2Species[,i]
	names(temp)<-row.names(mappingResult2Species)
	temp<-rev(sort(temp))
	if (length(temp)>10) {
		temp<-c(temp[1:9],sum(temp[10:length(temp)]))
		names(temp)[10]<-"Other"	
	}
	pie(temp,col=rainbow(10),cex=0.5)
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

png(paste0(resultFile,".top",sigNum,".png"),width=1400,height=700,res=300)
ggplot(mappingResult2SpeciesSelectedForFigure,aes(x=Sample,y=Reads,fill=Species))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()




