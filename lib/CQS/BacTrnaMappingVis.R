resultFile<-commandArgs()[7]
trnaCountTableFile<-commandArgs()[8]

trnaCountTable<-read.delim(trnaCountTableFile,header=T,row.names=1)

trnaCountTableSplit<-NULL
temp1<-strsplit(row.names(trnaCountTable),";")
temp1<-sapply(temp1,length)
trnaCountTableSplit<-trnaCountTable/temp1

nameUngroup<-strsplit(row.names(trnaCountTableSplit),";")
nameUngroup3<-sapply(nameUngroup,function(x) sapply(strsplit(x,"\\|"), function(x) x[3]))
nameUngroup7<-sapply(nameUngroup,function(x) sapply(strsplit(x,"\\|"), function(x) x[7]))

for ( i in 1:ncol(trnaCountTableSplit)) {
	png(paste0(resultFile,"_",colnames(trnaCountTableSplit)[i],".tRNA.png"),height=2000,width=4000,res=300)
	par(mfrow=c(1,2))
	par(mar=c(2,6,2,6))
	
	temp<-nameUngroup3
	nameToCount<-rep(0,length(unique(unlist(temp))))
	names(nameToCount)<-unique(unique(unlist(temp)))
	trnaCountTableSplitOneSample<-trnaCountTableSplit[,i]
	for (x in 1:length(temp)) {
		nameToCount[temp[[x]]]<-nameToCount[temp[[x]]]+trnaCountTableSplitOneSample[x]
	}
	nameToCountForFigure<-na.omit(rev(sort(nameToCount)))
	nameToCountForFigure<-c(nameToCountForFigure[1:9],Other=sum(nameToCountForFigure[-(1:9)]))
	pie(nameToCountForFigure,col=rainbow(10),main=paste0("Mapped Reads: ",as.integer(sum(nameToCountForFigure))))
	
	temp<-nameUngroup7
	nameToCount<-rep(0,length(unique(unlist(temp))))
	names(nameToCount)<-unique(unique(unlist(temp)))
	trnaCountTableSplitOneSample<-trnaCountTableSplit[,i]
	for (x in 1:length(temp)) {
		nameToCount[temp[[x]]]<-nameToCount[temp[[x]]]+trnaCountTableSplitOneSample[x]
	}
	nameToCountForFigure<-na.omit(rev(sort(nameToCount)))
	nameToCountForFigure<-c(nameToCountForFigure[1:9],Other=sum(nameToCountForFigure[-(1:9)]))
	pie(nameToCountForFigure,col=rainbow(10),main=paste0("Mapped Reads: ",as.integer(sum(nameToCountForFigure))))
	
	dev.off()
}

