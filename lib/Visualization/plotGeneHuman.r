library(ggplot2)
library(data.table)
library(dplyr)
library(Homo.sapiens)
library(ggbio)
library(GenomicRanges)
library(biovizBase)
library(cowplot)

data(genesymbol, package = "biovizBase")

if(! exists("inputFile")){
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
  
  locus=unique(dataForPlot$Locus)
  
  fileNumber=length(unique(dataForPlot$FileName))
  xlim=c(floor(min(dataForPlot$Position) / 1000) * 1000, ceiling(max(dataForPlot$Position) / 1000) * 1000)

  if(selectedFeature != locus){
    wh <- genesymbol[c(selectedFeature)]
    wh <- range(wh, ignore.strand = FALSE)
    g2<-autoplot(Homo.sapiens, which = wh) + 
      xlim(xlim) + 
      theme_classic() + 
      theme(axis.line=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin = margin(0, 1, 0, 0, "cm"))
  }

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
    pdf(outputFile, height=max(6, length(unique(dataForPlot$FileName))), width=8, onefile=TRUE)
    g <- ggplot(dataForPlot, aes_string("Position", column)) +
      ylab(columns[column]) +
      xlim(xlim) +
      geom_line()

    if(featureName != locus){
      g <- g + ggtitle(paste0(featureName, " (", locus, ")"))
    }

    if(hasCNV){
      g <- g + geom_area(aes(fill=CNV), alpha=.5)+
        scale_fill_manual(values=cols)
    }else{
      g <- g + geom_area()
    }

    g <- g +
      facet_grid(cols="FileName") +
      theme_bw() +
      theme(strip.background=element_blank(),
            axis.title.x=element_blank())

    if(selectedFeature != locus){
      g<-plot_grid(g, g2@ggplot, ncol = 1, align = 'v', axis = "lr", rel_heights=c(fileNumber, 1) )
    }
    
    print(g)
    dev.off()
  }
}

