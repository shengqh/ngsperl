resultFile<-outFile
trnaCountTableFile<-parFile1
groupFileList<-parSampleFile1
#resultFile<-commandArgs()[7]
#trnaCountTableFile<-commandArgs()[8]

library(reshape)
library(ggplot2)
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
groupPie<-function(x,maxCategory=10) {
	nameToCountForFigure<-na.omit(rev(sort(x)))
	if (length(nameToCountForFigure)>maxCategory) {
		nameToCountForFigure<-c(nameToCountForFigure[1:(maxCategory-1)],Other=sum(nameToCountForFigure[-(1:(maxCategory-1))]))
	}
	pie(nameToCountForFigure,col=rainbow(length(nameToCountForFigure)),main=paste0("Mapped Reads: ",as.integer(sum(nameToCountForFigure))))
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


trnaCountTable<-read.delim(trnaCountTableFile,header=T,row.names=1)
sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)

trnaCountTableExpand<-expandCountTableByName(trnaCountTable)

nameSub<-strsplit(row.names(trnaCountTableExpand),"_tRNA-|.trna\\d+-")
nameSubSpecies<-sapply(nameSub,function(x) x[1])
nameSubSpecies12<-sapply(strsplit(nameSubSpecies,"_"),function(x) paste0(x[1:2],collapse="_"))
temp<-sapply(nameSub,function(x) x[2])
temp<-gsub("-\\d+-\\d+$","",temp)
nameSubtRNA<-paste0(substr(temp,0,3),"-",substr(temp,nchar(temp)-2,nchar(temp)),sep="")
nameSubtRNA1<-substr(temp,0,3)

trnaCountTableExpandBySpecies<-aggregateCountTable(trnaCountTableExpand,nameSubSpecies)
trnaCountTableExpandBySpecies12<-aggregateCountTable(trnaCountTableExpand,nameSubSpecies12)
trnaCountTableExpandByRNA<-aggregateCountTable(trnaCountTableExpand,nameSubtRNA)
trnaCountTableExpandByRNA1<-aggregateCountTable(trnaCountTableExpand,nameSubtRNA1)


for ( i in 1:ncol(trnaCountTable)) {
	png(paste0(resultFile,"_",colnames(trnaCountTable)[i],".tRNA.png"),height=3000,width=3000,res=300)
	par(mfrow=c(2,2))
	par(mar=c(2,2,2,2))
	
	temp<-as.matrix(trnaCountTableExpandBySpecies)
	groupPie(temp[,i])
	
	temp<-as.matrix(trnaCountTableExpandBySpecies12)
	groupPie(temp[,i])
	
	temp<-as.matrix(trnaCountTableExpandByRNA)
	groupPie(temp[,i])
	
	temp<-as.matrix(trnaCountTableExpandByRNA1)
	groupPie(temp[,i])
	dev.off()
}

png(paste0(resultFile,".top.png"),width=3000,height=1500,res=300)
groupBarplot(trnaCountTableExpandByRNA1,groupName="tRNA")
dev.off()
