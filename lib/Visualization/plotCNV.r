library(ggplot2)
library(data.table)

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

rawTable<-fread(file=inputFile, sep="\t", header=T)

cols <- c("DEL" = "blue", "DUP" = "red", "NoIndel" = "gray")
column<-"PositionCount"
for (column in c("Percentage", "PositionCount")){
  selectedFeature<-unique(rawTable$Feature)[1]
  for (selectedFeature in unique(rawTable$Feature)) {
    cat(paste0(selectedFeature, "\n"))
    dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]

    outputFile = paste0(outputPrefix, ".", selectedFeature, ".", column, ".pdf")
    if (file.exists(outputFile)){
      next
    }
    
    pdf(outputFile, height=max(10, length(unique(dataForPlot$File))), width=10, onefile=TRUE)
    g <- ggplot(dataForPlot, aes_string("Position", column)) + 
      ggtitle(selectedFeature) +
      geom_line()
    
    if(("CNV" %in% colnames(dataForPlot)) & (length(unique(dataForPlot$CNV)) > 1)){
      dataForPlot$CNV[dataForPlot$CNV=="FALSE"] = "DEL"
      dataForPlot$CNV[dataForPlot$CNV=="NOREAD"] = "DEL"
      dataForPlot$CNV[dataForPlot$CNV==""] = "NoIndel"
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
