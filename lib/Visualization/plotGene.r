library(ggplot2)
library(data.table)
library(dplyr)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    setwd("/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_10_CNVGenesPlot/result")
    inputFile<-"linton_exomeseq_3321.position.txt.1000"
    outputPrefix<-'linton_exomeseq_3321.position'
    sizeFactorFile<-"/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_08_SizeFactor/result/linton_exomeseq_3321.txt.sizefactor"
  }else{
    inputFile<-args[1]
    outputPrefix<-args[2]
    sizeFactorFile<-args[3]
  }
}

cat(inputFile, "\n")
rawTable<-fread(file=inputFile, sep="\t", header=T)

readSizeFactors<-function(fileName){
  sizeFactorsAll<-read.table(sizeFactorFile, sep="\t", header=T)
  sizeFactorsAll<-sizeFactorsAll[!is.na(sizeFactorsAll$SizeFactor),]
  sfMedian<-sizeFactorsAll %>% group_by(Sample) %>% summarise(Median=median(SizeFactor))

  result<-sfMedian$Median
  names(result)<-sfMedian$Sample
  return(result)
}

sizeFactors<-readSizeFactors(sizeFactorFile)

cols <- c("DEL" = "blue", "DUP" = "red", "NoIndel" = "gray")
selectedFeature<-unique(rawTable$Feature)[1]
for (selectedFeature in unique(rawTable$Feature)) {
  cat(paste0(selectedFeature, "\n"))

  dataForPlot<-rawTable[which(rawTable$Feature==selectedFeature),]
  dataForPlot$NormalizedPositionCount<-round(dataForPlot$PositionCount * sizeFactors[dataForPlot$File])
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
      ylab(columns[column]) +
      xlab(unique(dataForPlot$Locus))+ 
      geom_line() 

    if(featureName != unique(dataForPlot$Locus)){
      g <- g + ggtitle(featureName)
    }

    if(hasCNV){
      g <- g + geom_area(aes(fill=CNV), alpha=.5)+
        scale_fill_manual(values=cols)
    }else{
      g <- g + geom_area(fill = "black")
    }

    g <- g +
      facet_grid(rows=FileName~.) +
      theme_bw() +
      theme(strip.background=element_blank())

    print(g)
    dev.off()
  }
}
