library(ggplot2)

args <- commandArgs(TRUE)
if (length(args) >= 2) {
  positionFile<-args[1]
  outputPrefix<-args[2]
}else{ #debug
  positionFile<-"/nobackup/vickers_lab/projects/20250321_Pua_12444_smallRNA/data_visualization/host_genome_yRNA_position_vis/result/Pua_12444.yRNAhost.position"
  outputPrefix<-"/nobackup/vickers_lab/projects/20250321_Pua_12444_smallRNA/data_visualization/host_genome_yRNA_position_vis/result/Pua_12444.position"
}

rawTable<-read.delim(positionFile,header=T,as.is=T,check.names=F)

selectedFeature=rawTable$Feature[1]
for (selectedFeature in unique(rawTable$Feature)) {
  dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]
  
  selectedFeature_fn<-gsub('\\?', '',selectedFeature)
  selectedFeature_fn<-gsub(':', '_',selectedFeature_fn)
  selectedFeature_fn<-gsub('\\|', '_',selectedFeature_fn)

  height = ceiling(length(unique(dataForPlot$File)) * 0.15 + 0.5)
  width = 8

  g = ggplot(dataForPlot, aes(Position, File)) + 
    geom_tile(aes(fill = Percentage)) + 
    scale_fill_gradientn(colours=rev(heat.colors(100))) + 
    ggtitle(selectedFeature) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=12, face="bold"),
          axis.text.y = element_text(face="bold"),
          plot.title = element_text(hjust = 0.5, size=14, face="bold"))
  ggsave(paste0(selectedFeature_fn,".Percentage.pdf"), g, height=height, width=width)
  ggsave(paste0(selectedFeature_fn,".Percentage.png"), g, height=height, width=width, units="in", dpi=300, bg="white")

  g = ggplot(dataForPlot, aes(Position, File)) +
    geom_tile(aes(fill = PositionCount)) + 
    scale_fill_gradientn(colours=rev(heat.colors(100))) + 
    ggtitle(selectedFeature) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=12, face="bold"),
          axis.text.y = element_text(face="bold"),
          plot.title = element_text(hjust = 0.5, size=14, face="bold"))
  ggsave(paste0(selectedFeature_fn,".PositionCount.pdf"), g, height=height, width=width)
  ggsave(paste0(selectedFeature_fn,".PositionCount.png"), g, height=height, width=width, units="in", dpi=300, bg="white")
}

