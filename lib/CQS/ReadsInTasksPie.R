options(bitmapType='cairo')

#############################
#Count reads in all NonParallel Tasks, Pie Chart;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
hostFile<-parFile1
nonhostFile<-parFile2
groupFileList<-parSampleFile1
groupVisLayoutFileList<-parSampleFile2

if(!exists("visLayoutAlphabet")){
  visLayoutAlphabet<-FALSE
}

facetColCount=getFacetColCount(groupFileList)

hostTable<-read.delim(hostFile,header=T,row.names=1,comment.char = '#')
hostTable<-hostTable[c("FeatureReads", "GenomeReads", "TooShortReads"),,drop=F]
rownames(hostTable)<-c('Host Small RNA','Mapped to Host Genome','Too Short for Mapping')

nonhostTable<-read.delim(nonhostFile,header=T,row.names=1,comment.char = '#')
if("UnmappedReads" %in% rownames(nonhostTable)){
  nonhostTable<-nonhostTable[c("FeatureReads", "UnmappedReads"),,drop=F]
}else{
  nonhostTable<-nonhostTable[c("FeatureReads", "UnannotatedReads"),,drop=F]
}
rownames(nonhostTable)<-c('Mapped to Non-Host','Unmapped')

tableForPieChart<-rbind(hostTable, nonhostTable)
write.csv(tableForPieChart,paste0(resultFile,".NonParallel.TaskReads.csv"))

ggpieToFile(tableForPieChart,fileName=paste0(resultFile,".NonParallel.TaskReads.Piechart.png"),maxCategory=NA,textSize=textSize,reOrder=FALSE,facetColCount=facetColCount)
#Group Pie chart
ggpieGroupToFile(tableForPieChart,fileName=paste0(resultFile,".NonParallel.TaskReads.Group.Piechart.png"),maxCategory=NA,reOrder=FALSE,
    groupFileList=groupFileList,
    outFileName=paste0(resultFile,".NonParallel.TaskReads.PercentGroups.csv"),textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList,
    visLayoutAlphabet=visLayoutAlphabet)
