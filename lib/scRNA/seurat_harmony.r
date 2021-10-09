rm(list=ls()) 
outFile='scRNA_6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/scratch/cqs/shengq2/annet_kirabo_projects/20210809_6383_scRNA_human/seurat_rawdata/result/scRNA_6383.rawobj.rds'
parFile2=''
parFile3=''


setwd('C:/projects/scratch/cqs/shengq2/annet_kirabo_projects/20210809_6383_scRNA_human/seurat_harmony/result')

source("scRNA_func.r")


library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(cowplot)
library(scales)
library(stringr)
library(harmony)
library(patchwork)
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

prefix<-outFile

has_batch_file<-file.exists(parSampleFile2)
npcs<-as.numeric(myoptions$pca_dims)
pca_dims<-1:npcs

rawobj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(rawobj, myoptions, prefix)
rawobj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

objs<-SplitObject(object = rawobj, split.by = "sample")
rm(rawobj)

obj=do_harmony(objs, by_sctransform, npcs, parSampleFile2)
reduction="harmony"
rm(objs)

for (reduct in c("pca", "harmony")){
  png(paste0(outFile, ".elbowplot.", reduct, ".png"), width=1500, height=1200, res=300)
  p<-ElbowPlot(obj, ndims = 20, reduction = reduct)
  print(p)
  dev.off()
}

cat("run_umap ... ")
obj <- RunUMAP(object = obj, reduction=reduction, dims=pca_dims, verbose = FALSE)

finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
saveRDS(finalList, file=finalListFile)

output_integraion_dimplot(obj, outFile, has_batch_file)
