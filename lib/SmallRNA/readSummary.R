rm(list=ls()) 
outFile='nextflex'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/preprocessing/fastqc_post_trim_summary/result/nextflex.countInFastQcVis.Result.Reads.csv'
parFile2='/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/data_visualization/reads_in_tasks/result/nextflex.NonParallel.TaskReads.csv'
parFile3='/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/host_genome/bowtie1_genome_1mm_NTA_smallRNA_category/result/nextflex.Category.Table.csv'


setwd('/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/data_visualization/read_summary/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
options(bitmapType='cairo')

myoptions = read_file_map(parSampleFile1, do_unlist = FALSE)

resultPrefix<-outFile
allreadsFile<-parFile1
readCategoryFile<-parFile2
smallRNACategoryFile<-parFile3

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(grid)

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

tableBarplot<-function(dat,textSize=13,colorNames="Set1",barwidth=0.5, legend_space_cm=0.5, ylab="Reads") {
  datForFigure<-getDatForFigure(dat)
  if(length(colorNames)==1){
    colors<-makeColors(length(unique(datForFigure[,"Category"])),colorNames)
  }else{
    colors<-colorNames
  }
  p<-ggplot(datForFigure) +
    geom_bar(aes(x=Sample,y=Reads,fill=Category), stat="identity", width=barwidth) +
    scale_fill_manual(values=colors) +
    theme_classic() +
    facet_grid(cols="Type", scale="free_y") +
    theme(legend.position = "top")+
    guides(fill = guide_legend(nrow = 1)) +
    ylab(ylab) +
    theme(axis.text.x=element_text(face="bold", size=textSize, angle=90,hjust=1,vjust=0.5),
          axis.text.y=element_text(size=textSize),
          legend.text=element_text(face="bold", size=textSize),
          legend.spacing.x = unit(legend_space_cm, 'cm'),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face="bold", size=textSize),
          legend.title= element_blank(),
          strip.text = element_text(face="bold", size=textSize),
          strip.background = element_blank())
  
  return(p)
}

truncNames=function(x,ncharMax=20) {
	if(max(nchar(x))>ncharMax) {
		x=make.unique(substr(x,0,ncharMax))
	}
	return(x)
}

hasReadCategory <- readCategoryFile != ''

allreads<-read.csv(allreadsFile, row.names=1, check.names=F)
allreads_name_map=c("Reads for Mapping"="Reads for Mapping", "Removed by Removing Sequence"="Removed by Sequence", "Removed by Trimming"="Removed by Trimming")
rownames(allreads) <- allreads_name_map[rownames(allreads)]
colnames(allreads)=truncNames(colnames(allreads)) #make names shorter for figure if they are too long

smallRNACategory<-read.csv(smallRNACategoryFile, row.names=1, check.names=F)
smallRNACategory<-smallRNACategory[grepl("ERV|RNA", rownames(smallRNACategory)),]

colnames(smallRNACategory)=truncNames(colnames(smallRNACategory))

g1<-tableBarplot(allreads, legend_space_cm=1, ylab="Reads\n")

smallRNA_colors=get_host_smallRNA_colors(myoptions$color_theme)

if(myoptions$color_theme == 'vickers') {
  name_map=c("ERV"="ERV", "lncRNA"="lincDR", "miRNA"="miRNA", "misc_RNA"="osRNA", "mt_tRNA"="mt_tDR", "rRNA"="rDR", "snoRNA"="snoDR", "snRNA"="snDR", "tRNA"="tDR", "yRNA"="yDR")
  rownames(smallRNACategory) <- name_map[rownames(smallRNACategory)]
  names(smallRNA_colors) <- name_map[names(smallRNA_colors)]
}

g3<-tableBarplot(smallRNACategory, colorNames=smallRNA_colors, ylab="Assigned reads\n")

p1 <- ggplot_gtable(ggplot_build(g1))
dev.off()
p3 <- ggplot_gtable(ggplot_build(g3))
dev.off()

if(hasReadCategory) {
  readCategory<-read.csv(readCategoryFile, row.names=1)
  if(myoptions$color_theme == 'vickers') {
    host_rows=c("Host Small RNA", "Mapped to Host Genome")
    other_rows=c("Mapped to Non-Host", "Unmapped", "Too Short for Mapping")
    readCategory <- rbind(colSums(readCategory[host_rows,]),
                          readCategory[other_rows,])
    rownames(readCategory) <- c("Host", other_rows)
  }
  readCategory_name_map=c("Host"="Host", "Host Small RNA"="Host smallRNA", "Mapped to Host Genome"="Host Genome", "Too Short for Mapping"="TooShort", "Mapped to Non-Host"="Non-host", "Unmapped"="Unknown")
  rownames(readCategory) <- readCategory_name_map[rownames(readCategory)]
  colnames(readCategory)=truncNames(colnames(readCategory))
}else{
  mappedReads<-apply(smallRNACategory, 2, sum)
  unmappedReads<-allreads["Reads for Mapping",] - mappedReads
  readCategory<-data.frame("Host smallRNA"=mappedReads, "Unknown" = unlist(unmappedReads), check.names=F)
  readCategory<-data.frame(t(readCategory))
}
category_colors=get_category_reads_colors(myoptions$color_theme)
g2<-tableBarplot(readCategory, colorNames=category_colors, ylab="Assigned reads\n")
p2 <- ggplot_gtable(ggplot_build(g2))
dev.off()

# maxWidth = unit.pmax(p1$widths[2:3], p3$widths[2:3])
# maxWidth = unit.pmax(maxWidth, p2$widths[2:3])
# p2$widths[2:3] <- maxWidth
# p1$widths[2:3] <- maxWidth
# p3$widths[2:3] <- maxWidth

maxHeight = 1600

width<-max(4000, 60 * ncol(allreads))
png(file=paste0(resultPrefix, ".png"), width=width, height=3 * maxHeight, res=300)
gg1<-ggarrange(p1, p2, p3, ncol = 1, nrow=3, labels=c("A", "B", "C"))
print(gg1)
dev.off()

