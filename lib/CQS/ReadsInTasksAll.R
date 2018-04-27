options(bitmapType='cairo')

#############################
#Count reads in all Tasks, including Table tasks and summrization in ReadsInTaskPie;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
readsInTaskCountFile<-parFile1
readsInTaskPieSummaryFile<-parFile2
totalCountFile<-parFile3

readsInTaskCount<-read.csv(readsInTaskCountFile,header=T,as.is=T,row.names=1)
readsInTaskPieSummary<-read.csv(readsInTaskPieSummaryFile,header=T,as.is=T,row.names=1)

readsInTaskAll<-rbind(readsInTaskCount,readsInTaskPieSummary[c("Mapped to Host Genome","Too Short for Mapping","Unmapped"),])


tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,"All.TaskReads.Barplot.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Category")
tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,"All.TaskReads.PerMillion.Barplot.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Category")
tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,"All.TaskReads.Barplot2.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"))
tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,"All.TaskReads.PerMillion.Barplot2.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"))

