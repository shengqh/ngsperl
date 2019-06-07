library(ggplot2)
library(data.table)

args <- commandArgs(TRUE)
if(length(args) == 0){
  setwd("/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_9_CNVGenesPlot/result")
  inputFile<-"linton_exomeseq_3321.position.CCL3L3.txt"
  outputPrefix<-'linton_exomeseq_3321.position'
  sizeFactorFile<-"/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/background/linton_exomeseq_3321.excluded.bed.sizefactor"
  #inputFile<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position.txt'
  #outputPrefix<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position'
}else{
  inputFile<-args[1]
  outputPrefix<-args[2]
  sizeFactorFile<-args[3]
}

rawTable<-fread(file=inputFile, sep="\t", header=T)

#sizeFactors<-read.table(sizeFactorFile, sep="\t", header=T)

cols <- c("DEL" = "blue", "DUP" = "red", "NoIndel" = "gray")
column<-"PositionCount"
selectedFeature<-unique(rawTable$Feature)[1]
for (selectedFeature in unique(rawTable$Feature)) {
  cat(paste0(selectedFeature, "\n"))
  dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]

  hasCNV<-(("CNV" %in% colnames(dataForPlot)) & (length(unique(dataForPlot$CNV)) > 1))
  
  if(hasCNV){
    dataForPlot$CNV[dataForPlot$CNV=="FALSE"] = "DEL"
    dataForPlot$CNV[dataForPlot$CNV=="NOREAD"] = "DEL"
    dataForPlot$CNV[dataForPlot$CNV==""] = "NoIndel"
  }
  
  for (column in c("Percentage", "PositionCount")){
    outputFile = paste0(outputPrefix, ".", selectedFeature, ".", column, ".pdf")
    #if (file.exists(outputFile)){
    #  next
    #}
    pdf(outputFile, height=max(10, length(unique(dataForPlot$File))), width=10, onefile=TRUE)
    g <- ggplot(dataForPlot, aes_string("Position", column)) + 
      ggtitle(selectedFeature) +
      geom_line()
    
    if(hasCNV){
      g <- g + geom_area(aes(fill=CNV), alpha=.5)+
        scale_fill_manual(values=cols)
    }
    
    g <- g +
      facet_grid(File~., scales = "free_y") +
      theme_bw() +
      theme(strip.background=element_blank())
    
    print(g)
    dev.off()
  }
}
