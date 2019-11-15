if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    setwd("/scratch/jbrown_lab/shengq2/projects/20191101_chipseq_4074_mouse/macs1callpeak_homer_annotation/result/ILEUM_FXR/")
    inputFile<-"ILEUM_FXR_peaks.name.bed.annotation.stats"
    outputPrefix<-'ILEUM_FXR_peaks.name.bed.annotation.stats'
  }else{
    inputFile<-args[1]
    outputPrefix<-args[2]
  }
}

cat(inputFile, "\n")

options(bitmapType='cairo')
library(ggplot2)
library(ggpubr)

stats<-read.table(inputFile, sep="\t", header=T, stringsAsFactors=F)
detailAnnotationIndex = which(stats$Annotation == "Annotation")
stats<-stats[c(1:(detailAnnotationIndex - 1)),c(1:2)]
colnames(stats)<-c("Category", "PeakNumber")
stats$PeakNumber<-as.numeric(stats$PeakNumber)
stats<-stats[stats$PeakNumber > 0,]
stats<-stats[order(stats$PeakNumber, decreasing=T),]
stats$Category<-factor(stats$Category, levels=(stats$Category))
stats$Percentage<-stats$PeakNumber * 100 / sum(stats$PeakNumber)

write.csv(stats, file=paste0(outputPrefix, ".csv"), row.names=F, quote=F)

legend.font<-element_text(size=14)
pdf(file=paste0(outputPrefix, ".pdf"), width=6, height=5)
g<-ggplot(stats, aes(x="", y=PeakNumber, fill=Category)) + geom_col() + coord_polar(theta = "y", direction = -1) + theme_void() +
	theme(legend.title=legend.font, legend.text=legend.font, legend.margin=margin(c(0,20,0,0)))
print(g)
dev.off()
