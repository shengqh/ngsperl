#############################
#Count reads in all table Tasks;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
countFilesList<-parSampleFile1

countFiles<-read.delim(countFilesList,header=F,as.is=T)

resultTable<-NULL
for (countFile in countFiles[,1]) {
	if (grepl(".csv$",countFile)) {
		countTable<-read.csv(countFile,header=T,row.names=1)
	} else {
		countTable<-read.delim(countFile,header=T,row.names=1)
	}
	resultTable<-rbind(resultTable,colSums(countTable))
}


row.names(resultTable)<-basename(countFiles[,1])

tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.Barplot.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,
		fill=NA,facet="Category")
