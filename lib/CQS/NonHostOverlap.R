options(bitmapType='cairo')
library(plyr)

#############################
#Vis for all Non Host Reads: Group 1, 2, 4; tRNA; two rRNA Categry;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
readFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3
totalCountFile<-parFile3

facetColCount=getFacetColCount(groupFileList)

categoriesNames<-c("Microbiome","Environment","Fungus","tRNA","rRNA")

readsMappingNames<-list()
readsMappingTable<-NULL
readFiles<-read.delim(readFileList,header=F,as.is=T)
for (i in 1:nrow(readFiles)) {
  cat("reading ", readFiles[i,1], "...\n")
  temp<-read.delim(readFiles[i,1],header=T,row.names=1,as.is=T)
  readsMappingNames[[i]]<-row.names(temp)
  readsMappingTable<-rbind.fill(readsMappingTable,temp)
}
readsMappingTable[is.na(readsMappingTable)]<-0
######################################
#Reads Overlap: Reads were found in how many categories? For reads only in one category, which one?
######################################
readsMappingNamesTable<-table(unlist(readsMappingNames))
dataForPlot<-NULL
#maxReadCategory<-max(readsMappingNamesTable)
#for (i in 2:maxReadCategory) {
#	temp<-colSums(readsMappingTable[names(readsMappingNamesTable)[which(readsMappingNamesTable==i)],])
#	dataForPlot<-rbind(dataForPlot,temp)
#}
#row.names(dataForPlot)<-paste0("Reads Mapped to ",2:maxReadCategory," Categories")

readsInOneCategory<-names(readsMappingNamesTable)[which(readsMappingNamesTable==1)]
for (i in 1:length(readsMappingNames)) {
	temp<-intersect(readsInOneCategory,readsMappingNames[[i]])
	temp<-colSums(readsMappingTable[temp,])
	dataForPlot<-rbind(dataForPlot,temp)
}
temp<-colSums(readsMappingTable[names(readsMappingNamesTable)[which(readsMappingNamesTable>1)],])
dataForPlot<-rbind(dataForPlot,temp)
#row.names(dataForPlot)[(nrow(dataForPlot)-length(categoriesNames)+1):nrow(dataForPlot)]<-categoriesNames
row.names(dataForPlot)<-c(paste0(categoriesNames," Only"),"Mapped to more than one Category")

write.csv(dataForPlot,paste0(resultFile,".Overlap.csv"))

#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".Piechart.png"),maxCategory=maxCategory,textSize=textSize, facetColCount=facetColCount)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


######################################
#Reads Overlap: Venn for 5 categories
######################################
if (length(readsMappingNames)>5) {
	dataForPlot<-readsMappingNames[1:5]
} else {
	dataForPlot<-readsMappingNames
}
names(dataForPlot)<-categoriesNames[1:5]
colors<-makeColors(length(categoriesNames))
vennCex=1.2
for (i in 1:ncol(readsMappingTable)) {
	reads2count<-readsMappingTable[,i]
	names(reads2count)<-row.names(readsMappingTable)
	png(paste0(resultFile,".",colnames(readsMappingTable)[i],".venn.png"),res=300,height=2000,width=2000)
	venn.diagram1(dataForPlot,count=reads2count,cex=vennCex,cat.cex=vennCex,fill=colors,alpha=0.7,margin=0.2,cat.dist=c(0.2,0.3,0.2,0.2,0.3))
	dev.off()
}



###################################################
#Group1 and Group2 reads mapping table overlap
###################################################
temp1<-intersect(readsMappingNames[[1]],readsMappingNames[[2]])
temp2<-setdiff(readsMappingNames[[1]],temp1)
temp3<-setdiff(readsMappingNames[[2]],temp1)
resultOut<-rbind(readsMappingTable[temp1,],readsMappingTable[temp2,],readsMappingTable[temp3,])
resultOut<-cbind(Category=c(rep("BothCategories",length(temp1)),rep("MicrobiomeOnly",length(temp2)),rep("EnvironmentOnly",length(temp3))),resultOut)
write.csv(resultOut,paste0(resultFile,".MicrobiomeVsEnvironment.reads.csv"))

temp1<-colSums(readsMappingTable[temp1,])
temp2<-colSums(readsMappingTable[temp2,])
temp3<-colSums(readsMappingTable[temp3,])
dataForPlot<-rbind(BothCategories=temp1,MicrobiomeOnly=temp2,EnvironmentOnly=temp3)

#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Piechart.png"),maxCategory=maxCategory,textSize=textSize,facetColCount=facetColCount)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


