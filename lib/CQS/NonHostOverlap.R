#############################
#Vis for all Non Host Reads: Group 1, 2, 4; tRNA; two rRNA Categry;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
readFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3
totalCountFile<-parFile3

readsMappingNames<-list()
readsMappingTable<-NULL
readFiles<-read.delim(readFileList,header=F,as.is=T)
for (i in 1:nrow(readFiles)) {
	temp<-read.delim(readFiles[i,1],header=T,row.names=1,as.is=T)
	readsMappingNames[[i]]<-row.names(temp)
	readsMappingTable<-rbind(readsMappingTable,temp)
}

#Reads Overlap: Reads were found in how many categories?
readsMappingNamesTable<-table(unlist(readsMappingNames))
dataForPlot<-NULL
maxReadCategory<-max(readsMappingNamesTable)
for (i in 1:maxReadCategory) {
	temp<-colSums(readsMappingTable[names(readsMappingNamesTable)[which(readsMappingNamesTable==i)],])
	dataForPlot<-rbind(dataForPlot,temp)
}
row.names(dataForPlot)<-paste0("Reads Mapped to ",1:maxReadCategory," Categories")

#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)

###################################################
#Group1 and Group2 reads mapping table overlap
###################################################
temp1<-intersect(readsMappingNames[[1]],readsMappingNames[[2]])
temp2<-setdiff(readsMappingNames[[1]],temp1)
temp3<-setdiff(readsMappingNames[[2]],temp1)

temp1<-colSums(readsMappingTable[temp1,])
temp2<-colSums(readsMappingTable[temp2,])
temp3<-colSums(readsMappingTable[temp3,])
dataForPlot<-rbind(BothCategories=temp1,MicrobiomeOnly=temp2,EnvironmentOnly=temp3)

#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


