library(ggplot2)

args <- commandArgs(TRUE)
if(length(args) == 0){
  inputFile<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position.txt'
  outputPrefix<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position'
}else{
  inputFile<-args[1]
  outputPrefix<-args[2]
}

rawTable<-read.delim(inputFile,header=T,as.is=T,stringsAsFactor=F)

for (column in c("Percentage", "PositionCount")){
  outputFile = paste0(outputPrefix, ".", column, ".pdf")
  pdf(outputFile, height=10, width=10, onefile=TRUE)
  for (selectedFeature in unique(rawTable$Feature)) {
    dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]
  
    g <- ggplot(dataForPlot, aes(Position, File)) + 
         geom_tile(aes_string(fill = column)) + 
         scale_fill_gradientn(colours=rev(heat.colors(100))) + 
         ggtitle(selectedFeature)
    print(g)
  }
  dev.off()
}

