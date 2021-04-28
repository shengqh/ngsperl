options(bitmapType='cairo')
library(plyr)
rbindFillWithName<-function(x,y,...) {
	uniqueNamesY<-setdiff(row.names(y),row.names(x))
	if (length(uniqueNamesY)==0) {
		return(x)
	} else {
		allNames<-c(row.names(x),uniqueNamesY)
		result<-rbind.fill(x,y[uniqueNamesY,,drop=FALSE],...)
		row.names(result)<-allNames
		return(result)
	}
}

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

if(! exists("categoriesNames")){
  categoriesNames<-c("Microbiome","Environment","Fungus","tRNA","rRNA")
}

readsMappingNames<-list()
readsMappingTable<-NULL
readFiles<-read.delim(readFileList,header=F,as.is=T)
for (i in 1:nrow(readFiles)) {
  cat("reading ", readFiles[i,1], "...\n")
  temp<-read.delim(readFiles[i,1],header=T,row.names=1,as.is=T,check.names=F)
  readsMappingNames[[i]]<-row.names(temp)
  readsMappingTable<-rbindFillWithName(readsMappingTable,temp)
}
readsMappingTable[is.na(readsMappingTable)]<-0

######################################
#Reads Overlap: Reads were found in how many categories? For reads only in one category, which one?
######################################

cat("Preparing overlap reads now ","...\n")
allReadsNames<-unlist(readsMappingNames)
readsMappingNamesTable<-table(allReadsNames)
dataForPlot<-NULL
#maxReadCategory<-max(readsMappingNamesTable)
#for (i in 2:maxReadCategory) {
# temp<-colSums(readsMappingTable[names(readsMappingNamesTable)[which(readsMappingNamesTable==i)],])
# dataForPlot<-rbind(dataForPlot,temp)
#}
#row.names(dataForPlot)<-paste0("Reads Mapped to ",2:maxReadCategory," Categories")

readsInOneCategory<-names(readsMappingNamesTable)[which(readsMappingNamesTable==1)]
i<-1
for (i in 1:length(readsMappingNames)) {
  readsOne<-intersect(readsInOneCategory,readsMappingNames[[i]])
  readsOneTable<-readsMappingTable[intersect(allReadsNames,readsOne),]
  readsOneColSums<-colSums(readsOneTable)
  dataForPlot<-rbind(dataForPlot,readsOneColSums)
}
readsInMoreCategory<-names(readsMappingNamesTable)[which(readsMappingNamesTable>1)]
readsMoreTable<-readsMappingTable[intersect(allReadsNames,readsInMoreCategory),]
readsMoreColSums<-colSums(readsMoreTable)
dataForPlot<-rbind(dataForPlot,readsMoreColSums)
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
cat("Making Venn diagram for overlap reads now ","...\n")
if (length(readsMappingNames)>5) {
   warning("More than 5 categories. Only first 5 categories will be used in overlap Venn diagram.")
   dataForPlot<-readsMappingNames[1:5]
 } else {
   dataForPlot<-readsMappingNames
 }
 names(dataForPlot)<-categoriesNames[1:5]
# colors<-makeColors(length(categoriesNames))
 colors<-makeColors(5)
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
common12<-intersect(readsMappingNames[[1]],readsMappingNames[[2]])
cat1only<-setdiff(readsMappingNames[[1]],common12)
cat2only<-setdiff(readsMappingNames[[2]],common12)

reads12MappingTable<-readsMappingTable[c(1:(length(readsMappingNames[[1]]) + length(readsMappingNames[[2]]))),]
reads12MappingNames = c(readsMappingNames[[1]], readsMappingNames[[2]])

common12table<-reads12MappingTable[reads12MappingNames %in% common12,]
cat1onlytable<-reads12MappingTable[reads12MappingNames %in% cat1only,]
cat2onlytable<-reads12MappingTable[reads12MappingNames %in% cat2only,]

resultOut<-rbind(common12table, cat1onlytable, cat2onlytable)
resultOut<-data.frame(Read=c(reads12MappingNames[reads12MappingNames %in% common12], reads12MappingNames[reads12MappingNames %in% cat1only], reads12MappingNames[reads12MappingNames %in% cat2only]), 
                 Category=c(rep("Both_Microbiome",length(common12)),
                            rep("Both_Environment",length(common12)),
                            rep("MicrobiomeOnly",nrow(cat1onlytable)),
                            rep("EnvironmentOnly",nrow(cat2onlytable))),resultOut,stringsAsFactors=FALSE)
resultOut<-resultOut[-which(resultOut$Category=="Both_Environment"),]
resultOut$Category<-gsub("Both_Microbiome","Both",resultOut$Category)
write.csv(resultOut,paste0(resultFile,".MicrobiomeVsEnvironment.reads.csv"),row.names=F)

common12colsums<-colSums(common12table)
cat1onlycolsums<-colSums(cat1onlytable)
cat2onlycolsums<-colSums(cat2onlytable)
dataForPlot<-rbind(BothCategories=common12colsums,MicrobiomeOnly=cat1onlycolsums,EnvironmentOnly=cat2onlycolsums)

#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Piechart.png"),maxCategory=maxCategory,textSize=textSize,facetColCount=facetColCount)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Group.Piechart.png"),groupFileList=groupFileList,
    outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


###################################################
#Reads Overlap: Venn for Group1, Group2 and Group4
###################################################
cat("Making Venn diagram for non host genome overlap reads now ","...\n")
dataForPlot<-readsMappingNames[1:3]
names(dataForPlot)<-categoriesNames[1:3]
colors<-makeColors(3)
vennCex=1.2
for (i in 1:ncol(readsMappingTable)) {
	reads2count<-readsMappingTable[,i]
	names(reads2count)<-row.names(readsMappingTable)
	png(paste0(resultFile,".",colnames(readsMappingTable)[i],".NonHostGenome.venn.png"),res=300,height=2000,width=2000)
	venn.diagram1(dataForPlot,count=reads2count,cex=vennCex,cat.cex=vennCex,fill=colors,alpha=0.7,margin=0.2,cat.dist=c(0.2,0.2,0.2))
	dev.off()
}

###################################################
#Reads Overlap: Venn for tRNA, and rRNA
###################################################
cat("Making Venn diagram for non host library overlap reads now ","...\n")
dataForPlot<-readsMappingNames[4:5]
names(dataForPlot)<-categoriesNames[4:5]
colors<-makeColors(3)[1:2]
vennCex=1.2
for (i in 1:ncol(readsMappingTable)) {
	reads2count<-readsMappingTable[,i]
	names(reads2count)<-row.names(readsMappingTable)
	png(paste0(resultFile,".",colnames(readsMappingTable)[i],".NonHostLibrary.venn.png"),res=300,height=2000,width=2000)
	venn.diagram1(dataForPlot,count=reads2count,cex=vennCex,cat.cex=vennCex,fill=colors,alpha=0.7,margin=0.2,cat.dist=c(0.2,0.2))
	dev.off()
}


