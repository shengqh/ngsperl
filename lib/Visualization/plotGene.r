library(ggplot2)
library(data.table)

args <- commandArgs(TRUE)
if(length(args) == 0){
  setwd("/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_10_CNVGenesPlot/result")
  inputFile<-"linton_exomeseq_3321.position.txt.1000"
  outputPrefix<-'linton_exomeseq_3321.position'
  sizeFactorFile<-"/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_08_SizeFactor/result/linton_exomeseq_3321.txt.sizefactor"
  #inputFile<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position.txt'
  #outputPrefix<-'/scratch/cqs/shengq2/vickers/20190504_smallRNA_as_chipseq_GCF_000005845.2_ASM584v2/plotPeak/result/Control.position'
}else{
  inputFile<-args[1]
  outputPrefix<-args[2]
  sizeFactorFile<-args[3]
}

cat(inputFile)
rawTable<-fread(file=inputFile, sep="\t", header=T)

sizeFactors<-read.table(sizeFactorFile, sep="\t", header=T)
rownames(sizeFactors)<-paste0(sizeFactors$Sample, ":", sizeFactors$Chromosome)

cols <- c("DEL" = "blue", "DUP" = "red", "NoIndel" = "gray")
selectedFeature<-unique(rawTable$Feature)[1]
for (selectedFeature in unique(rawTable$Feature)) {
  cat(paste0(selectedFeature, "\n"))

  dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]
  dataForPlot$Key<-paste0(dataForPlot$File, ":", dataForPlot$Chromosome)
  dataForPlot$NormalizedPositionCount<-round(dataForPlot$PositionCount * sizeFactors[dataForPlot$Key, "SizeFactor"])
  dataForPlot$FileName<-gsub("_", " ", dataForPlot$File)

  hasCNV<-(("CNV" %in% colnames(dataForPlot)) & (length(unique(dataForPlot$CNV)) > 1))
  
  if(hasCNV){
    dataForPlot$CNV[dataForPlot$CNV=="FALSE"] = "DEL"
    dataForPlot$CNV[dataForPlot$CNV=="NOREAD"] = "DEL"
    dataForPlot$CNV[dataForPlot$CNV==""] = "NoIndel"
  }
  
  if(grepl("_\\d+_\\d+$", selectedFeature)){
    parts=unlist(strsplit(selectedFeature, "_"))
    featureName=paste0(paste0(parts[1:(length(parts)-2)]), ":", parts[length(parts)-1], "-", parts[length(parts)])
  }else{
    featureName = selectedFeature
  }

  columns<-c("Read count", "Percentage", "Normalized read count")
  names(columns)<-c("PositionCount", "Percentage", "NormalizedPositionCount")
  column<-"NormalizedPositionCount"
  for (column in names(columns)){
    outputFile = paste0(outputPrefix, ".", selectedFeature, ".", column, ".pdf")
    #if (file.exists(outputFile)){
    #  next
    #}
    pdf(outputFile, height=max(4, length(unique(dataForPlot$FileName))), width=6, onefile=TRUE)
    g <- ggplot(dataForPlot, aes_string("Position", column)) + 
      ggtitle(featureName) +
      ylab(columns[column]) +
      geom_line()
    
    if(hasCNV){
      g <- g + geom_area(aes(fill=CNV), alpha=.5)+
        scale_fill_manual(values=cols)
    }
    
    g <- g +
      facet_grid(FileName~.) +
      theme_bw() +
      theme(strip.background=element_blank())
    
    print(g)
    dev.off()
  }
}
