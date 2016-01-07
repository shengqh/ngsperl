library("reshape2")
library("ggplot2")

#args<-commandArgs(trailingOnly=TRUE)
#inputFile<-args[1]
#outputFile<-args[2]

setwd("H:/shengquanhu/projects/JonathanBrown/20151215_janathan_atacseq_fat/macs2bdgdiff_replicates_nomodel_depth/result/SQ_Visc_CHOW")
inputFile<-"SQ_Visc_CHOW_c3.0_cond1.bed.depth"
outputFile<-"SQ_Visc_CHOW_c3.0_cond1.bed.depth.pdf"

data<-read.table(inputFile, sep="\t", header=T, stringsAsFactors = F)
files<-unique(data$File)

pdf(outputFile, onefile = T)
for(x in files){
  cat(x, "\n")
  curdata<-data[data$File==x,]
  mdata<-melt(curdata, id=c("Chr", "Position", "File"))
  colnames(mdata)<-c("Chr", "Position", "File", "Sample", "Depth")
  g<-ggplot(mdata, aes(x=Position, y=Depth, colour=Sample)) + 
    geom_point() + 
    stat_smooth() + 
    xlab(unique(data$chr)) + 
    ggtitle(x) +
    facet_grid( . ~ Sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  print(g)
}
dev.off()
