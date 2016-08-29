#############################
#Vis for all Non Host Reads: Group 1, 2, 4; tRNA; two rRNA Categry;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
readFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3
totalCountFile<-parFile3

readsMappingResult<-list()
for (i in 1:nrow(readFileList)) {
	readsMappingResult[[i]]<-read.delim(readFileList[i,1],header=T,row.names=1)
}
readsMappingResultTable<-sapply(readsMappingResult,function(x) x)

#Group1 and Group2 reads mapping table overlap
temp<-table(c(row.names(readTable1),row.names(readTable2)))
readsOverlapTime<-list()
for (i in 1:max(temp)) {
	readsOverlapTime[[i]]<-names(temp)[which(temp==i)]
}


overlapCount<-colSums(readTable1[readsOverlap,])
Genome1OnlyCount<-colSums(readTable1[setdiff(row.names(readTable1),readsOverlap),])
Genome2OnlyCount<-colSums(readTable2[setdiff(row.names(readTable2),readsOverlap),])

dataForPlot<-rbind(MicrobiomeOnly=Genome1OnlyCount,EnvironmentOnly=Genome2OnlyCount,InBoth=overlapCount)

#Pie chart for all samples
resultFile<-"Test"
totalCountFile<-""
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dat=dataForPlot,fileName=paste0(resultFile,".Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dat=mappingResult2Species,fileName=paste0(resultFile,".Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)
mappingResult2Species<-countTableToSpecies(dat=mappingResult,databaseLogFile=databaseLogFile,outFileName=paste0(resultFile,".Species.csv"))

#Pie chart for all samples
ggpieToFile(mappingResult2Species,fileName=paste0(resultFile,".Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dat=mappingResult2Species,fileName=paste0(resultFile,".Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dat=mappingResult2Species,fileName=paste0(resultFile,".Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)



