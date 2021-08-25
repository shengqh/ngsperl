options(bitmapType='cairo')

#############################
#Count reads in all table Tasks;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
countFilesList<-parSampleFile1
totalCountFile<-parFile3

countFiles<-read.delim(countFilesList,header=F,as.is=T)

resultTable<-NULL
taskFileWithReads<-NULL
for (countFile in countFiles[,1]) {
	if (grepl(".csv$",countFile)) {
		countTable<-read.csv(countFile,header=T,row.names=1,check.names=F)
	} else {
		countTable<-read.delim(countFile,header=T,row.names=1,check.names=F)
	}
	if (nrow(countTable)==0) {
		next;
	}
	colClass<-sapply(countTable, class)
	countTableNum<-countTable[,which(colClass=="numeric" | colClass=="integer")]
	
	if (is.null(resultTable)) {
		resultTable<-rbind(colSums(countTableNum))
	} else {
		resultTable<-rbind(resultTable,colSums(countTableNum)[colnames(resultTable)])
	}
	taskFileWithReads<-c(taskFileWithReads,countFile)
}

#change names
row.names(resultTable)<-gsub("_pm_.+","",basename(taskFileWithReads))
row.names(resultTable)<-gsub(".+.miRNA.count$","Host_genome_miRNA",row.names(resultTable))
row.names(resultTable)<-gsub(".+.tRNA.count$","Host_genome_tRNA",row.names(resultTable))
row.names(resultTable)<-gsub(".+.yRNA.count$","Host_genome_yRNA",row.names(resultTable))
row.names(resultTable)<-gsub(".+.snRNA.count$","Host_genome_snRNA",row.names(resultTable))
row.names(resultTable)<-gsub(".+.snoRNA.count$","Host_genome_snoRNA",row.names(resultTable))
row.names(resultTable)<-gsub(".+.rRNA.count$","Host_genome_rRNA",row.names(resultTable))
row.names(resultTable)<-gsub(".+.other.count$","Host_genome_other_smallRNA",row.names(resultTable))
write.csv(resultTable,paste0(resultFile,".TaskReads.csv"))

tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.Barplot.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Category",proportionBar=FALSE)
tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.PerMillion.Barplot.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Category",proportionBar=FALSE)
tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.Barplot2.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"),proportionBar=FALSE)
tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.PerMillion.Barplot2.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"),proportionBar=FALSE)

