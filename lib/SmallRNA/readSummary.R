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

tableBarplot<-function(dat,textSize=4,ylab="Reads",colorNames="Set1",barwidth=0.5, percent=FALSE) {
  dat[is.na(dat)]<-0
  if(percent){
    dat<-prop.table(as.matrix(dat), 2)
    ylab<-paste0("Percentage of ", ylab)
  }
  datForFigure<-melt(as.matrix(dat))
  colnames(datForFigure)<-c("Category","Sample","Reads")
  colors<-makeColors(length(unique(datForFigure[,"Category"])),colorNames)
  p<-ggplot(datForFigure) +
    geom_bar(aes(x=Sample,y=Reads,fill=Category), stat="identity", width=barwidth) +
    scale_fill_manual(values=colors) +
    theme_classic() +
    theme(legend.position = "top")+
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text = element_text(size=textSize),legend.text=element_text(size=textSize),
          axis.title = element_text(size=textSize),legend.title= element_text(size=textSize))+
    ylab(ylab)
  
  return(p)
}

allreads<-read.csv(allreadsFile, row.names=1)
readCategory<-read.csv(readCategoryFile, row.names=1)
smallRNACategory<-read.csv(smallRNACategoryFile, row.names=1)
smallRNACategory<-smallRNACategory[grepl("RNA", rownames(smallRNACategory)),]

for (perc in c(TRUE, FALSE)){
  g1<-tableBarplot(allreads, percent=perc)
  g2<-tableBarplot(readCategory, percent=perc)
  g3<-tableBarplot(smallRNACategory, percent=perc)

  p1 <- ggplot_gtable(ggplot_build(g1))
  p2 <- ggplot_gtable(ggplot_build(g2))
  p3 <- ggplot_gtable(ggplot_build(g3))

  maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3])
  p1$widths[2:3] <- maxWidth
  p2$widths[2:3] <- maxWidth
  p3$widths[2:3] <- maxWidth

  label1<-c("A", "B", "C")
  
  width<-max(1500, 20 * ncol(allreads))
  png(file=paste0(resultPrefix, ifelse(perc, ".perc.png", ".count.png")), width=width, height=2000, res=300)
  gg1<-ggarrange(p1, p2, p3, ncol = 1, nrow=3, labels=label1)
  print(gg1)
  dev.off()
}