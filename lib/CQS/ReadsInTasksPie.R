#############################
#Count reads in all NonParallel Tasks, Pie Chart;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
countFilesList<-parSampleFile1
#totalCountFile<-parFile3

countFiles<-read.delim(countFilesList,header=F,as.is=T)

taskName<-sapply(strsplit(countFiles[,1],"\\/"),function(x) {resultFolderInd<-grep("result",x);return(x[resultFolderInd-1])})
totalCountAll<-NULL
for (countFile in countFiles[,1]) {
	if (grepl(".csv$",countFile)) {
		countTable<-read.csv(countFile,header=T,row.names=1,as.is=T)
	} else {
		countTable<-read.delim(countFile,header=T,row.names=1,as.is=T)
	}
	countInd<-which(colnames(countTable) %in% c("EstimateCount","Count"))
	if (length(countInd)>0) {
		totalCount<-sum(countTable[,countInd[1]])
	} else {
		stop(paste0("Can't find colnames matched to Count data in file ",countFile)) 
	}
	totalCountAll<-c(totalCountAll,totalCount)
}
resultTable<-data.frame(Task=taskName,Count=totalCountAll,Sample=countFiles[,2])
resultTable<-acast(resultTable,Task~Sample,value.var="Count")
write.csv(resultTable,paste0(resultFile,".NonParallel.TaskReads.csv"))

HostSmallRnaMapped<-resultTable["bowtie1_genome_1mm_NTA_smallRNA_count",]
HostGenomeMappedReads<-resultTable["identical",]-resultTable["bowtie1_genome_unmapped_reads",]-resultTable["bowtie1_genome_1mm_NTA_smallRNA_count",]
NonHostMapped<-resultTable["bowtie1_genome_unmapped_reads",]-resultTable["bowtie1_unmapped_reads",]
UnMapped<-resultTable["bowtie1_unmapped_reads",]

tableForPieChart<-rbind(HostSmallRnaMapped,HostGenomeMappedReads,NonHostMapped,UnMapped)
ggpieToFile(tableForPieChart,fileName=paste0(resultFile,".NonParallel.TaskReads.Piechart.png"),maxCategory=NA,textSize=textSize,reOrder=FALSE)




