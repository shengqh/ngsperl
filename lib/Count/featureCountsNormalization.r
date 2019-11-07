library(ggplot2)
library(data.table)
library(dplyr)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    setwd("/scratch/cqs/shengq2/jennifer/20190906_lindsay_exomeseq_3772_hg38/bwa_refine_nosoftclip_featureCounts_table/result")
    inputFile<-"lindsay_exomeseq_3772.count"
    outputPrefix<-'lindsay_exomeseq_3772.count'
    sizeFactorFile<-"/scratch/cqs/shengq2/jennifer/20190906_lindsay_exomeseq_3772_hg38/bwa_refine_nosoftclip_gatk4_CNV_Germline_08_SizeFactor/result/lindsay_exomeseq_3772.txt.sizefactor"
  }else{
    inputFile<-args[1]
    outputPrefix<-args[2]
    sizeFactorFile<-args[3]
  }
}

cat(inputFile, "\n")
rawTable<-read.table(file=inputFile, sep="\t", header=T, row.names=1)

readSizeFactors<-function(fileName){
  sizeFactorsAll<-read.table(sizeFactorFile, sep="\t", header=T)
  sizeFactorsAll<-sizeFactorsAll[!is.na(sizeFactorsAll$SizeFactor),]
  sfMedian<-sizeFactorsAll %>% group_by(Sample) %>% summarise(Median=median(SizeFactor))

  result<-sfMedian$Median
  names(result)<-sfMedian$Sample
  return(result)
}

sizeFactors<-readSizeFactors(sizeFactorFile)

normData<-rawTable
for(col in colnames(normData)){
  normData[,col] = normData[,col] * sizeFactors[col]
}
