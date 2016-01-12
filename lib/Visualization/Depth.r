library("reshape2")
library("ggplot2")

args<-commandArgs(trailingOnly=TRUE)
singlePdf<-(args[1] == "1")
inputFile<-args[2]
outputFile<-args[3]

if(is.null(outputFile)){
  outputFile = ""
}

setwd("H:/shengquanhu/projects/JonathanBrown/20151215_janathan_atacseq_fat/macs2bdgdiff_replicates_nomodel_depth/result/SQ_Visc_CHOW")
inputFile<-"SQ_Visc_CHOW_c3.0_cond1.bed.depth"
singlePdf<-FALSE
outputFile<-""

data<-read.table(inputFile, sep="\t", header=T, stringsAsFactors = F)
files<-unique(data$File)[1:5]

if(singlePdf){
  pdf(outputFile, onefile = T)
}

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
    facet_grid( Sample ~ .) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  if(singlePdf){
    print(g)
  }else{
    png(paste0(outputFile, x, ".png"), width=2000, height=max(3000, 300+500 * length(unique(curdata$Sample))), res=300)
    print(g)
    dev.off()
  }
}

if(singlePdf){
  dev.off()
}
