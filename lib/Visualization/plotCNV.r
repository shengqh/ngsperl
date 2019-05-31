library(ggplot2)

args <- commandArgs(TRUE)
if(length(args) == 0){
  inputFile<-"/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321.position.txt"
  outputPrefix<-'/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_8_PlotGeneCNV/result/linton_exomeseq_3321.position'
  #inputFile<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position.txt'
  #outputPrefix<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position'
}else{
  inputFile<-args[1]
  outputPrefix<-args[2]
}

rawTable<-read.delim(inputFile,header=T,as.is=T,stringsAsFactor=F)

for (column in c("Percentage", "PositionCount")){
  outputFile = paste0(outputPrefix, ".", column, ".pdf")
  pdf(outputFile, height=max(10, length(unique(rawTable$File))), width=10, onefile=TRUE)
  for (selectedFeature in unique(rawTable$Feature)) {
    dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]
    g <- ggplot(dataForPlot, aes_string("Position", column)) + 
      ggtitle(selectedFeature) +
      geom_line() +
      facet_grid(File~., scales = "free_y") +
      theme_bw() +
      theme(strip.background=element_blank())
    print(g)
  }
  dev.off()
}
