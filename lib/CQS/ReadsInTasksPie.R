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

if(!exists("uniqueGroupNames")){
  uniqueGroupNames<-NULL
}

facetColCount=getFacetColCount(groupFileList)

hostTable<-read.delim(hostFile,header=T,row.names=1,comment.char = '#',check.names=F)
hostTable<-hostTable[c("FeatureReads", "GenomeReads", "TooShortReads"),,drop=F]
rownames(hostTable)<-c('Host Small RNA','Mapped to Host Genome','Too Short for Mapping')

nonhostTable<-read.delim(nonhostFile,header=T,row.names=1,comment.char = '#',check.names=F)
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


if(file.exists(parFile3) & file.exists(groupFileList)){
  library(ggplot2)
  library(reshape2)
  microbial<-read.delim(parFile3, header=T, row.names=1,check.names=F)
  all<-tableForPieChart
  all<-data.frame(t(all))
  all$all<-rowSums(all)
  all$Human<-(all$Host.Small.RNA + all$Mapped.to.Host.Genome) / all$all * 100
  all$Microbial<-(microbial[rownames(all), "Count"]) / all$all * 100
  all$Other<-100 - all$Human - all$Microbial
  figData<-all[,c("Human", "Microbial", "Other")]
  if(!is.null(uniqueGroupNames)){
    group<-read.delim(groupFileList, stringsAsFactors = F, header=F,check.names=F)
    group<-group[group$V2 %in% uniqueGroupNames,]
  }else{
    group<-getSampleInGroup(groupFileList, rownames(microbial), useLeastGroups=TRUE, onlySamplesInGroup=TRUE)
  }
  rownames(group)<-group$V1
  figData$Group<-group[rownames(figData), "V2"]
  figData=figData[!is.na(figData$Group),]

  figData$Sample<-rownames(figData)
  figData<-figData[order(figData$Microbial, decreasing = T),]
  mFigData<-melt(figData, id.vars=c("Sample", "Group"))
  colnames(mFigData)<-c("Sample", "Group", "Category", "Percentage")
  mFigData$Sample<-factor(mFigData$Sample, levels=figData$Sample)
  mFigData$Category<-factor(mFigData$Category, levels=c("Microbial", "Host", "Other"))
  if(!is.null(uniqueGroupNames)){
    mFigData$Group<-factor(mFigData$Group, levels=uniqueGroupNames)
  }else{
    mFigData$Group<-factor(mFigData$Group)
  }
  pdf(paste0(outFile,".NonParallel.TaskReads.bar.pdf"), width=10, height=6)
  colors<-c("Microbial" = "chartreuse3", "Host" = "deepskyblue", "Other" = "gray")
  g<-ggplot(mFigData) + geom_bar(aes(y = Percentage, x = Sample, fill = Category), stat="identity") + facet_grid(~Group, scales = "free_x") +
    scale_fill_manual(values=colors) +
    ylab("Percentage of Reads") +
    ylim(0, 100) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background=element_blank())
  print(g)
  dev.off()
}
