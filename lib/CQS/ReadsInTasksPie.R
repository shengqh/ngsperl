rm(list=ls()) 
outFile='nextflex'
parSampleFile1='fileList1_pie.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_nextflex.mapped.count'
parFile2='/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/final_unmapped/final_unmapped_reads_summary/result/nextflex.count'
parFile3='/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/nonhost_genome/nonhost_genome_count/result/nextflex.nonhost_genome.tsv'
textSize=9;groupTextSize=10; 

setwd('/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/data_visualization/reads_in_tasks/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
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

readCategory2 <- rbind(colSums(tableForPieChart[1:2,]),
                       tableForPieChart[c(4,5,3),])
rownames(readCategory2) <- c("Host", "Non-host", "Unknown", "TooShort")

doGetDatForFigure<-function(dat,percent){
  dat[is.na(dat)]<-0
  if(percent){
    dat<-prop.table(as.matrix(dat), 2)
  }
  result<-reshape2::melt(as.matrix(dat))
  colnames(result)<-c("Category","Sample","Reads")
  result$Type=ifelse(percent, "Percentage", "Reads")
  return(result)
}

getDatForFigure<-function(dat){
  d1<-doGetDatForFigure(dat,percent=T)
  d2<-doGetDatForFigure(dat,percent=F)
  result<-rbind(d1,d2)
  return(result)
}

tableBarplot<-function(dat,textSize=13,colorNames="Set1",barwidth=0.5) {
  datForFigure<-getDatForFigure(dat)
  colors<-colorNames#makeColors(length(unique(datForFigure[,"Category"])),colorNames)
  p<-ggplot(datForFigure) +
    geom_bar(aes(x=Sample,y=Reads,fill=Category), stat="identity", width=barwidth) +
    scale_fill_manual(values=colors) +
    theme_classic() +
    facet_grid(Type~., scale="free_y") +
    theme(legend.position = "top")+
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size=textSize),legend.text=element_text(size=textSize),
          axis.title = element_text(size=textSize),legend.title= element_text(size=textSize),
          strip.background = element_blank()) +
    theme(strip.text.y = element_blank(),
          axis.text.y=element_text(face="bold", size=13, color = "black"),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, face="bold", size=13, color = "black"),
          axis.title.y = element_text(size=13, face="bold", colour = "black"),
          legend.text = element_text(face = "bold", size = 13),
          legend.title=element_blank()) +
    xlab("") + ylab("Assigned Reads")
  
  return(p)
}

g=tableBarplot(readCategory2, colorNames = c("#ad07e3", "#ff0066", "#107f80", "black"))
width<-max(2000,60*ncol(readCategory2) + 100) / 300
ggsave(paste0(outFile,".NonParallel.TaskReads.bar.2.png"), g, width=width, height=6, dpi=300, units="in", bg="white")

if(file.exists(parFile3) & file.exists(groupFileList)){
  microbial<-read.delim(parFile3, header=T, row.names=1,check.names=F)
  all<-tableForPieChart
  all<-data.frame(t(all))
  all$all<-rowSums(all)
  all$Host<-(all$Host.Small.RNA + all$Mapped.to.Host.Genome) / all$all * 100
  all$Microbial<-(microbial[rownames(all), "Count"]) / all$all * 100
  all$Other<-100 - all$Host - all$Microbial
  figData<-all[,c("Host", "Microbial", "Other")]
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
  mFigData<-reshape2::melt(figData, id.vars=c("Sample", "Group")) |>
    dplyr::rename(Category=variable, Percentage=value)
  mFigData$Sample<-factor(mFigData$Sample, levels=figData$Sample)
  mFigData$Category<-factor(mFigData$Category, levels=c("Microbial", "Host", "Other"))
  if(!is.null(uniqueGroupNames)){
    mFigData$Group<-factor(mFigData$Group, levels=uniqueGroupNames)
  }else{
    mFigData$Group<-factor(mFigData$Group)
  }
  colors<-c("Microbial" = "chartreuse3", "Host" = "deepskyblue", "Other" = "gray")
  g<-ggplot(mFigData) + 
    geom_bar(aes(y = Percentage, x = Sample, fill = Category), stat="identity") + 
    facet_grid(~Group, scales = "free_x") +
    scale_fill_manual(values=colors) +
    ylab("Proportion") +
    ylim(0, 100) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(face="bold", angle=90, hjust=1, vjust=0.5),
          axis.ticks.x=element_blank(),
          strip.background=element_blank())

  ggsave(paste0(outFile,".NonParallel.TaskReads.bar.png"), g, width=width, height=3, dpi=300, units="in", bg="white")
}

