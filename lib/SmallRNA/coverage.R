library(ggplot2)

args <- commandArgs(TRUE)
positionFile<-args[1]
outputPrefix<-args[2]

rawTable<-read.delim(positionFile,header=T,as.is=T,check.names=F)

for (selectedFeature in unique(rawTable$Feature)) {
  dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]
  
  selectedFeature<-gsub('\\?', '',selectedFeature)
  selectedFeature<-gsub(':', '_',selectedFeature)
  selectedFeature<-gsub('\\|', '_',selectedFeature)
  
  pdf(paste0(outputPrefix, selectedFeature,".Percentage.pdf"),height=10, width=10)
  print(ggplot(dataForPlot, aes(Position, File)) + geom_tile(aes(fill = Percentage))+scale_fill_gradientn(colours=rev(heat.colors(100))))
  dev.off()
  
  pdf(paste0(outputPrefix, selectedFeature,".PositionCount.pdf"),height=10, width=10)
  print(ggplot(dataForPlot, aes(Position, File)) +geom_tile(aes(fill = PositionCount))+scale_fill_gradientn(colours=rev(heat.colors(100))))
  dev.off()
}

