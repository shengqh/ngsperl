rm(list=ls()) 
outFile='obesity'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/cqs/ravi_shah_projects/20220221_obesity_smallRNA_hg38/preprocessing/fastqc_post_trim_summary/result/obesity.countInFastQcVis.Result.Reads.csv'
parFile2=''
parFile3='/scratch/cqs/ravi_shah_projects/20220221_obesity_smallRNA_hg38/host_genome/bowtie1_genome_1mm_NTA_smallRNA_category/result/obesity.Category.Table.csv'


setwd('/scratch/cqs/ravi_shah_projects/20220221_obesity_smallRNA_hg38/data_visualization/read_summary/result')

### Parameter setting end ###

options(bitmapType='cairo')

DEBUG<-FALSE

if (DEBUG){
  resultPrefix<-"/data/stein_lab/mjo_sRNA_data/20180424_michelle_smallRNA_2868_human_320to475/data_visualization/read_summary/result/P2868_treatment"
  allreadsFile<-"/data/stein_lab/mjo_sRNA_data/20180424_michelle_smallRNA_2868_human_320to475/preprocessing/fastqc_post_trim/result/P2868_treatment.countInFastQcVis.Result.Reads.csv"
  readCategoryFile<-"/data/stein_lab/mjo_sRNA_data/20180424_michelle_smallRNA_2868_human_320to475/data_visualization/reads_in_tasks/result/P2868_treatment.NonParallel.TaskReads.csv"
  smallRNACategoryFile<-"/data/stein_lab/mjo_sRNA_data/20180424_michelle_smallRNA_2868_human_320to475/host_genome/bowtie1_genome_1mm_NTA_smallRNA_category/result/P2868_treatment.Category.Table.csv"
}else{
  resultPrefix<-outFile
  allreadsFile<-parFile1
  readCategoryFile<-parFile2
  smallRNACategoryFile<-parFile3
}

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(grid)

makeColors<-function(n,colorNames="Set1") {
  maxN<-brewer.pal.info[colorNames,"maxcolors"]
  if (n<=maxN) {
    colors<-brewer.pal(n, colorNames)
  } else {
    colors<-colorRampPalette(brewer.pal(maxN, colorNames))(n)
  }
  return(colors)
}

doGetDatForFigure<-function(dat,percent){
  dat[is.na(dat)]<-0
  if(percent){
    dat<-prop.table(as.matrix(dat), 2)
  }
  result<-melt(as.matrix(dat))
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

tableBarplot<-function(dat,textSize=4,colorNames="Set1",barwidth=0.5) {
  datForFigure<-getDatForFigure(dat)
  colors<-makeColors(length(unique(datForFigure[,"Category"])),colorNames)
  p<-ggplot(datForFigure) +
    geom_bar(aes(x=Sample,y=Reads,fill=Category), stat="identity", width=barwidth) +
    scale_fill_manual(values=colors) +
    theme_classic() +
    facet_grid(Type~., scale="free_y") +
    theme(legend.position = "top")+
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text = element_text(size=textSize),legend.text=element_text(size=textSize),
          axis.title = element_text(size=textSize),legend.title= element_text(size=textSize),
          strip.background = element_blank()) +
          xlab("") + ylab("")
  
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
colnames(allreads)=truncNames(colnames(allreads)) #make names shorter for figure if they are too long

smallRNACategory<-read.csv(smallRNACategoryFile, row.names=1, check.names=F)
smallRNACategory<-smallRNACategory[grepl("RNA", rownames(smallRNACategory)),]
colnames(smallRNACategory)=truncNames(colnames(smallRNACategory))

g1<-tableBarplot(allreads)
g3<-tableBarplot(smallRNACategory)

p1 <- ggplot_gtable(ggplot_build(g1))
dev.off()
p3 <- ggplot_gtable(ggplot_build(g3))
dev.off()

maxWidth = unit.pmax(p1$widths[2:3], p3$widths[2:3])
maxHeight = 1000
if(hasReadCategory) {
  readCategory<-read.csv(readCategoryFile, row.names=1)
  colnames(readCategory)=truncNames(colnames(readCategory))
}else{
  mappedReads<-apply(smallRNACategory, 2, sum)
  unmappedReads<-allreads["Reads for Mapping",] - mappedReads
  readCategory<-data.frame("Host smallRNA"=mappedReads, "Unmapped" = unlist(unmappedReads), check.names=F)
  readCategory<-data.frame(t(readCategory))
}
g2<-tableBarplot(readCategory)
p2 <- ggplot_gtable(ggplot_build(g2))
dev.off()
maxWidth = unit.pmax(maxWidth, p2$widths[2:3])
p2$widths[2:3] <- maxWidth
maxHeight = 1600

p1$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth

width<-max(2000, 40 * ncol(allreads))
png(file=paste0(resultPrefix, ".png"), width=width, height=2 * maxHeight, res=300)
gg1<-ggarrange(p1, p2, p3, ncol = 1, nrow=3, labels=c("A", "B", "C"))
print(gg1)
dev.off()


