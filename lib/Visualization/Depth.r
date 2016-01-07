library("reshape2")
library("ggplot2")

args<-commandArgs(trailingOnly=TRUE)
inputFile<-args[1]
outputFile<-args[2]
#setwd("H:/shengquanhu/projects/JonathanBrown/20151215_janathan_atacseq_fat")
#inputFile<-"SQ_Visc_CHOW_common_1.depth"
#outputFile<-"SQ_Visc_CHOW_common_1.depth.png"

data<-read.table(inputFile, sep="\t", header=T)
mdata<-melt(data, id=c("chr", "position"))

colnames(mdata)<-c("Chr", "Position", "Sample", "Depth")

png(outputFile, width=max(2000, 200 + 500 * (ncol(data) - 2)), height=2000, res=300)
g<-ggplot(mdata, aes(x=Position, y=Depth, colour=Sample)) + 
  geom_point() + 
  stat_smooth() + 
  xlab(unique(data$chr)) + 
  facet_grid( . ~ Sample) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g)
dev.off()
