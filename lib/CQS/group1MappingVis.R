resultFile<-outFile
mappingResultFile<-parFile1
databaseLogFile<-parFile2

#resultFile<-commandArgs()[7]
#mappingResultFile<-commandArgs()[8]
#databaseLogFile<-commandArgs()[9]

library(ggplot2)
library(reshape)

mappingResult<-read.delim(mappingResultFile,header=T,row.names=1)
databaseLog<-read.delim(databaseLogFile,header=F,skip=1,as.is=T)
row.names(mappingResult)<-gsub("\\.\\d+","",row.names(mappingResult))

id2Species<-databaseLog[,1]
names(id2Species)<-gsub(".fna","",databaseLog[,2])

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
write.csv(mappingResult2Species,paste0(resultFile,".toGenome.csv"))

for (i in 1:ncol(mappingResult2Species)) {
	png(paste0(colnames(mappingResult2Species)[i],".group1.png"),width=3000,height=1500,res=300)
	par(mar=c(2,9,2,9))
	temp<-mappingResult2Species[,i]
	names(temp)<-row.names(mappingResult2Species)
	temp<-rev(sort(temp))
	temp<-c(temp[1:9],sum(temp[10:length(temp)]))
	names(temp)[10]<-"Other"	
	pie(temp,col=rainbow(10),cex=0.5)
	dev.off()
}

#find the most significant row and combine the result
sigNum<-5
temp<-apply(mappingResult2Species,2,function(x) rev(order(x))[1:sigNum])
mappingResult2SpeciesSelected<-mappingResult2Species[unique(as.vector(temp)),]
mappingResult2SpeciesSelected<-rbind(mappingResult2SpeciesSelected,Other=colSums(mappingResult2Species[-unique(as.vector(temp)),,drop=FALSE]))

mappingResult2SpeciesSelectedForFigure<-mappingResult2SpeciesSelected
mappingResult2SpeciesSelectedForFigure$Species<-row.names(mappingResult2SpeciesSelected)
mappingResult2SpeciesSelectedForFigure<-melt(mappingResult2SpeciesSelectedForFigure)
colnames(mappingResult2SpeciesSelectedForFigure)<-c("Species","Sample","Reads")

pdf(paste0(resultFile,".top",sigNum,".pdf"),width=14)
ggplot(mappingResult2SpeciesSelectedForFigure,aes(x=Sample,y=Reads,fill=Species))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()




