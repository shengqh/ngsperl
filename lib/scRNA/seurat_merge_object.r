rm(list=ls()) 
outFile='PH_combine'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/scratch/cqs/shengq2/paula_hurley_projects/20221115_scRNA_7467_benign_hg38/seurat_rawdata/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(digest)
library(patchwork)
library(sparseMatrixStats)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

sample_pattern= myoptions$sample_pattern

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
hemoglobinPattern=myoptions$hemoglobinPattern


#read raw count dat
filelist1<-read.table(parSampleFile1, header=F, stringsAsFactors = F)
rawobjs = list()
fidx=1
fileMap<-split(filelist1$V1, filelist1$V2)
fileTitle<-names(fileMap)[1]

for(fileTitle in names(fileMap)) {
  fileName  = fileMap[[fileTitle]]
  cat(fileTitle, "\n")
  sobj<-readRDS(fileName)
  if(is.list(sobj)){
    sobj<-sobj$rawobj
  }
  
  if(sample_pattern != ""){
    cells = colnames(sobj)[grepl(sample_pattern, sobj$orig.ident)]
    sobj<-subset(sobj, cells=cells)
  }

  if(!('percent.mt' %in% colnames(sobj@meta.data))){
    sobj<-PercentageFeatureSet(object=sobj, pattern=Mtpattern, col.name="percent.mt")
  }

  if(!('percent.ribo' %in% colnames(sobj@meta.data))){
    sobj<-PercentageFeatureSet(object=sobj, pattern=rRNApattern, col.name = "percent.ribo")
  }

  if(!('percent.hb' %in% colnames(sobj@meta.data))){
    sobj<-PercentageFeatureSet(object=sobj, pattern=hemoglobinPattern, col.name="percent.hb")
  }   

  if(!('sample' %in% colnames(sobj@meta.data))){
    sobj$sample=sobj$orig.ident
  }

  rawobjs[[fileTitle]] = sobj
  rm(sobj)
}

if(length(rawobjs) == 1){
  rawobj <- rawobjs[[1]]
}else{
  rawobj <- merge(rawobjs[[1]], y = unlist(rawobjs[2:length(rawobjs)]), project = "integrated")
}
rm(rawobjs)

#rawobj<-readRDS("PH_combine.rawobj.rds")
output_rawdata(rawobj, outFile, Mtpattern, rRNApattern, hemoglobinPattern)
writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
