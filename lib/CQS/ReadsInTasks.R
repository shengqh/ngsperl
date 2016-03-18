#############################
#Count reads in all table Tasks;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
countFilesList<-parSampleFile1
totalCountFile<-parFile3

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
write.csv(resultTable,paste0(resultFile,".TaskReads.csv"))

tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.Barplot.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,
		fill=NA,facet="Category")
tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.PerMillion.Barplot.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,
		fill=NA,facet="Category")
tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.Barplot2.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"))
tableBarplotToFile(dat=resultTable,fileName=paste0(resultFile,".TaskReads.PerMillion.Barplot2.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"))

