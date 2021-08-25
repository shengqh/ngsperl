options(bitmapType='cairo')

#############################
#Count reads in all Tasks, including Table tasks and summrization in ReadsInTaskPie;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
readsInTaskCountFile<-parFile1
readsInTaskPieSummaryFile<-parFile2
totalCountFile<-parFile3

readsInTaskCount<-read.csv(readsInTaskCountFile,header=T,as.is=T,row.names=1,check.names=F)
readsInTaskPieSummary<-read.csv(readsInTaskPieSummaryFile,header=T,as.is=T,row.names=1,check.names=F)

readsInTaskAll<-rbind(readsInTaskCount,readsInTaskPieSummary[c("Mapped to Host Genome","Too Short for Mapping","Unmapped"),])
write.csv(readsInTaskAll,paste0(resultFile,".All.TaskReads.csv"))


tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,".All.TaskReads.Barplot.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Category",proportionBar=FALSE)
tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,".All.TaskReads.PerMillion.Barplot.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Category",proportionBar=FALSE)
tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,".All.TaskReads.Barplot2.png"),
		totalCountFile="",maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"),proportionBar=FALSE)
tableBarplotToFile(dat=readsInTaskAll,fileName=paste0(resultFile,".All.TaskReads.PerMillion.Barplot2.png"),
		totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize,height=2500,
		fill=NA,facet="Sample",x="Category",y="Reads",varName=c("Category","Sample","Reads"),proportionBar=FALSE)

