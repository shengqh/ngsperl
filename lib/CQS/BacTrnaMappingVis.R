resultFile<-outFile
trnaCountTableFile<-parFile1
#resultFile<-commandArgs()[7]
#trnaCountTableFile<-commandArgs()[8]

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


trnaCountTable<-read.delim(trnaCountTableFile,header=T,row.names=1)

trnaCountTableExpand<-expandCountTableByName(trnaCountTable)

nameSub<-strsplit(row.names(trnaCountTableExpand),"-")
nameSubSpecies<-sapply(nameSub,function(x) gsub("_tRNA","",x[1]))
nameSubSpecies12<-sapply(strsplit(nameSubSpecies,"_"),function(x) paste0(x[1:2],collapse="_"))
nameSubtRNA<-sapply(nameSub,function(x) paste0(x[2:3],collapse="-"))
nameSubtRNA1<-sapply(nameSub,function(x) x[2])

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
